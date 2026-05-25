import csv
import time

# BURAYA ZİPTEN ÇIKARDIĞIN DOSYANIN TAM ADINI YAZ
INPUT_FASTA = "rna.fna" 
OUTPUT_CSV = "database.csv"

print(f"Dönüşüm başlatılıyor... Hedef dosya: {INPUT_FASTA}")
start_time = time.time()

try:
    with open(INPUT_FASTA, "r", encoding="utf-8") as f_in, open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f_out:
        writer = csv.writer(f_out)
        # Sütun başlıklarımızı (Motorun beklediği başlıklar) yazıyoruz
        writer.writerow(["Gene_ID", "Sequence", "Description", "Biotype"])

        current_id = ""
        current_desc = ""
        current_seq = []
        count = 0

        for line in f_in:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith(">"):
                # Eski geni csv'ye kaydet
                if current_id:
                    # DNA'yı RNA formatına (T -> U) çevirerek kaydediyoruz (Opsiyonel ama daha sağlıklı)
                    full_seq = "".join(current_seq).replace("T", "U")
                    writer.writerow([current_id, full_seq, current_desc, "transcript"])
                    count += 1

                # Yeni genin bilgilerini parçala
                # Örnek NCBI Header: >NM_000014.6 Homo sapiens alpha-2-macroglobulin...
                parts = line[1:].split(" ", 1)
                current_id = parts[0]
                current_desc = parts[1] if len(parts) > 1 else "Unknown Target"
                current_seq = []
            else:
                # Sadece harfleri listeye ekle
                current_seq.append(line)

        # En son kalan geni de kaydet
        if current_id:
            full_seq = "".join(current_seq).replace("T", "U")
            writer.writerow([current_id, full_seq, current_desc, "transcript"])
            count += 1

    end_time = time.time()
    print(f"\nMUHTEŞEM! İşlem tamamlandı.")
    print(f"Toplam {count} adet gen dizilimi başarıyla database.csv dosyasına aktarıldı.")
    print(f"Geçen Süre: {round(end_time - start_time, 2)} saniye.")

except FileNotFoundError:
    print(f"HATA: '{INPUT_FASTA}' adlı dosya backend klasöründe bulunamadı!")
    print("Lütfen zipten çıkan FASTA dosyasının adını koda doğru yazdığına emin ol.")