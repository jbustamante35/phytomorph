 CMD = ['xargs -P 10 -I % -d ''\n'' -a fileName.txt fastOrderintoSets.sh %'];
 
 locate -d /mnt/scratch4/spaldingdata.db -d /mnt/scratch4/tetra.db -d /mnt/scratch4/spaldingimages.db -d /mnt/scratch4/snapper.db --regex (png|tif|jpg|bmp)$ > out.txt
 
 locate -d /mnt/scratch4/spaldingdata.db -d /mnt/scratch4/tetra.db -d /mnt/scratch4/spaldingimages.db -d /mnt/scratch4/snapper.db --regex ''' EXT '$'''
 
 locate -n 10 -r tif$ > out.txt
 
 
 
 xargs -P 12 -I % -d '\n' -a ./out.txt ./fastOrderintoSets.sh % ./pathNames.txt ./fileNames.txt

 sort -u pathNames.txt > sPathNames.txt
 nl -s '*' sPathNames.txt > IpathNames.txt
 
 
 ./superScanner.sh '1*/home/nate/' ./junkTest/out.txt
 ./fullScanner.sh ./junkTest/IpathNames.txt ./junkTest/out.txt
 
 
 locate -d /mnt/scratch4/spaldingdata.db -d /mnt/scratch4/tetra.db -d /mnt/scratch4/spaldingimages.db -d /mnt/scratch4/snapper.db --regex '(png|tif|jpg|bmp)$' | xargs -P 12 -I % -d '\n' ./fastOrderintoSets.sh % ./pathNames.txt ./fileNames.txt
 
 
 locate -d /mnt/scratch4/spaldingdata.db -d /mnt/scratch4/tetra.db -d /mnt/scratch4/spaldingimages.db -d /mnt/scratch4/snapper.db --regex '(png|tif|jpg|bmp)$' > totalFileList.txt
 xargs -a totalFileList.txt -n 100 -P 12 -I % -d '\n' ./fastOrderintoSets.sh % ./pathNames.txt ./fileNames.txt
 
 
 