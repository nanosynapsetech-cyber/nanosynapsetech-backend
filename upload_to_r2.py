#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# upload_to_r2.py
#
# Uploads kmer_output/ and seq_output/ to Cloudflare R2.
# Features:
#   - Resume support: skips files already in R2 (checks existing keys)
#   - Parallel uploads (configurable thread count)
#   - Retry on failure
#   - Progress bar
#
# Setup (run once):
#   Set the 4 variables below (R2_ACCOUNT_ID, R2_ACCESS_KEY, etc.)
#
# Usage:
#   python3 upload_to_r2.py

import os
import sys
import time
import json
import boto3
from botocore.config import Config
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock

# ─── CONFIG — FILL THESE IN ──────────────────────────────────────────────────
R2_ACCOUNT_ID  = os.environ.get("R2_ACCOUNT_ID",  "")   # Cloudflare Account ID
R2_ACCESS_KEY  = os.environ.get("R2_ACCESS_KEY",  "")   # R2 Access Key ID
R2_SECRET_KEY  = os.environ.get("R2_SECRET_KEY",  "")   # R2 Secret Access Key
R2_KMER_BUCKET = os.environ.get("R2_KMER_BUCKET", "nanosynapse-kmer")
R2_SEQ_BUCKET  = os.environ.get("R2_SEQ_BUCKET",  "nanosynapse-seq")

LOCAL_KMER_DIR = "kmer_output"
LOCAL_SEQ_DIR  = "seq_output"

MAX_WORKERS    = 8      # parallel upload threads (R2 drops conn if too high)
MAX_RETRIES    = 5      # retry on failure
UPLOAD_LOG     = "r2_upload_log.json"  # tracks uploaded files for resume

# ─── CHECKS ──────────────────────────────────────────────────────────────────
if not all([R2_ACCOUNT_ID, R2_ACCESS_KEY, R2_SECRET_KEY]):
    print("ERROR: Missing R2 credentials!")
    print()
    print("Set these environment variables before running:")
    print("  $env:R2_ACCOUNT_ID  = 'your-account-id'")
    print("  $env:R2_ACCESS_KEY  = 'your-access-key-id'")
    print("  $env:R2_SECRET_KEY  = 'your-secret-access-key'")
    print()
    print("Find them at: dash.cloudflare.com -> R2 -> Manage R2 API Tokens")
    sys.exit(1)

# ─── S3 CLIENT (Cloudflare R2 is S3-compatible) ──────────────────────────────
endpoint = f"https://{R2_ACCOUNT_ID}.r2.cloudflarestorage.com"

s3 = boto3.client(
    "s3",
    endpoint_url=endpoint,
    aws_access_key_id=R2_ACCESS_KEY,
    aws_secret_access_key=R2_SECRET_KEY,
    config=Config(
        retries={"max_attempts": MAX_RETRIES, "mode": "adaptive"},
        max_pool_connections=MAX_WORKERS + 5,
    ),
    region_name="auto",
)

# ─── HELPERS ─────────────────────────────────────────────────────────────────
_log_lock  = Lock()
_counter_lock = Lock()

def load_log() -> set:
    """Load set of already-uploaded file keys."""
    if os.path.exists(UPLOAD_LOG):
        with open(UPLOAD_LOG) as f:
            return set(json.load(f))
    return set()

def save_log(uploaded: set):
    with _log_lock:
        with open(UPLOAD_LOG, "w") as f:
            json.dump(sorted(uploaded), f)

def upload_file(bucket: str, local_path: str, key: str) -> bool:
    """Upload a single file. Returns True on success."""
    content_type = "application/json" if key.endswith(".json") else "text/plain"
    for attempt in range(MAX_RETRIES):
        try:
            s3.upload_file(
                local_path, bucket, key,
                ExtraArgs={"ContentType": content_type}
            )
            return True
        except Exception as e:
            if attempt < MAX_RETRIES - 1:
                time.sleep(2 ** attempt)
            else:
                print(f"\n  FAILED: {key} — {e}")
                return False
    return False

def upload_directory(local_dir: str, bucket: str, already_uploaded: set) -> tuple[int, int]:
    """
    Upload all files in local_dir to bucket.
    Returns (success_count, fail_count).
    """
    all_files = os.listdir(local_dir)
    to_upload = [f for f in all_files if f"{bucket}/{f}" not in already_uploaded]

    total   = len(all_files)
    skip    = total - len(to_upload)
    success = skip
    fail    = 0

    print(f"\n  Bucket   : {bucket}")
    print(f"  Total    : {total:,} files")
    print(f"  Skipping : {skip:,} already uploaded")
    print(f"  Uploading: {len(to_upload):,} new files")
    if not to_upload:
        print("  Nothing to upload!")
        return success, fail

    start = time.time()

    def _task(fname):
        local_path = os.path.join(local_dir, fname)
        key = fname
        ok = upload_file(bucket, local_path, key)
        return fname, ok

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(_task, f): f for f in to_upload}
        done = 0
        for future in as_completed(futures):
            fname, ok = future.result()
            done += 1
            with _counter_lock:
                if ok:
                    success += 1
                    already_uploaded.add(f"{bucket}/{fname}")
                else:
                    fail += 1

            if done % 500 == 0 or done == len(to_upload):
                elapsed = time.time() - start
                rate    = done / elapsed if elapsed > 0 else 0
                eta     = (len(to_upload) - done) / rate if rate > 0 else 0
                print(
                    f"    {done:>6,}/{len(to_upload):,}  "
                    f"({done/len(to_upload)*100:.1f}%)  "
                    f"{rate:.0f} files/s  "
                    f"ETA: {int(eta//60)}m{int(eta%60)}s",
                    end="\r", flush=True
                )
                # Save log every 500 files for resume support
                save_log(already_uploaded)

    print()  # newline after \r
    elapsed = time.time() - start
    print(f"  Done in {elapsed/60:.1f} min — success={success:,} fail={fail:,}")
    return success, fail


# ─── MAIN ────────────────────────────────────────────────────────────────────
print("=" * 55)
print("NanoSynapse — Cloudflare R2 Upload")
print("=" * 55)
print(f"Endpoint : {endpoint}")
print(f"Buckets  : {R2_KMER_BUCKET} + {R2_SEQ_BUCKET}")
print(f"Workers  : {MAX_WORKERS} parallel threads")

# Verify connection
print("\nVerifying R2 connection...")
try:
    s3.head_bucket(Bucket=R2_KMER_BUCKET)
    print(f"  {R2_KMER_BUCKET}: OK")
except Exception as e:
    print(f"  ERROR connecting to {R2_KMER_BUCKET}: {e}")
    print("  Make sure the bucket exists and credentials are correct.")
    sys.exit(1)

try:
    s3.head_bucket(Bucket=R2_SEQ_BUCKET)
    print(f"  {R2_SEQ_BUCKET}: OK")
except Exception as e:
    print(f"  ERROR connecting to {R2_SEQ_BUCKET}: {e}")
    sys.exit(1)

# Load resume log
uploaded = load_log()
print(f"\nResume log: {len(uploaded):,} files already uploaded")

total_start = time.time()

# Upload kmer_output/
print("\n[1/2] Uploading kmer_output/ ...")
kmer_ok, kmer_fail = upload_directory(LOCAL_KMER_DIR, R2_KMER_BUCKET, uploaded)
save_log(uploaded)

# Upload seq_output/
print("\n[2/2] Uploading seq_output/ ...")
seq_ok, seq_fail = upload_directory(LOCAL_SEQ_DIR, R2_SEQ_BUCKET, uploaded)
save_log(uploaded)

# Final summary
total_elapsed = time.time() - total_start
print(f"\n{'=' * 55}")
print("Upload Complete!")
print(f"  kmer_output: {kmer_ok:,} ok | {kmer_fail:,} failed")
print(f"  seq_output : {seq_ok:,} ok  | {seq_fail:,} failed")
print(f"  Total time : {total_elapsed/60:.1f} minutes")

if kmer_fail == 0 and seq_fail == 0:
    print("\nAll files uploaded successfully!")
    print("  Next: Configure main_turso.py to query R2")
else:
    print(f"\n{kmer_fail + seq_fail} files failed — re-run to retry (resume supported)")
