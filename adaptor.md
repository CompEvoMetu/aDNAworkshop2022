### Adaptör ne işe yarar?

-Dna fragmentlerini çoğaltmak işleminin gerçekleşeceği katı faza bağlamak \
-index ile bu fragmentlerin kimliğini oluşturmak 

$AdapterRemoval --file1 $1 --file2 $2 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNATCTCGTATGCCGTCTTCTGCTTG --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT --qualitybase 33 --gzip --qualitymax 60 --trimns --collapse --minalignmentlength 11 --threads ${Cores} --basename $3/$4 --settings $4.settings


--collapse--> Eşleştirilmiş uç modunda, iki çift çakışırsa, ikisini birleştirerek ve kalite puanlarını yeniden hesaplayarak iki okumayı tek bir okuma şeklinde birleştirir.\
--minalignmentlength--> Okumalar bire birleştirlmeden önce, eşleştirilmiş uç okumaları birleşetirirken veya tek uç modunda tam şablon dizilerini tanımlamaya çalışırken mate 1 ve mate 2 arasındaki minimum örtüşme (Default:11)

kırpılmış ve birleştirlmemiş okuma çiftlerini içeren fastq dosyaları \
*.pair1.truncated \
*.pair2.truncated \
*.singleton.truncated.gz birleştirilmrde atılan okulmar \
*.collapsed.gz birleştirilmiş okumalar \
*.collapsed.truncated.gz --trimns seçeneğiyle kırpılmış birleştirilmiş okumalar
