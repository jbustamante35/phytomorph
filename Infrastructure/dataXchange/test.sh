


#!/bin/bash
    

#<keyShare>
generateKey() {
 sudo rngd -r /dev/urandom
 echo -e "Key-Type: RSA\nKey-Length: 4096\nSubkey-Type: RSA\nSubkey-Length: 4096\nName-Real:tina.miller\nName-Comment:siteID_\nName-Email:this@works.com" > /tmp/junk.key
 gpg2 --batch --gen-key /tmp/junk.key
 rm /tmp/junk.key
 }
exportKey(){
 gpg2 --export "tina.miller" > /tmp/public.key
 }
sharePublicKey(){
 curl https://s3dev.chtc.wisc.edu/phytomorphservice/ -F bucket=phytomorphservice -F policy=eyJleHBpcmF0aW9uIjoiMjAyMC0wNS0wN1QyMjoxMzowNS4xMDNaIiwiY29uZGl0aW9ucyI6W1siZXEiLCIkYnVja2V0IiwicGh5dG9tb3JwaHNlcnZpY2UiXSxbInN0YXJ0cy13aXRoIiwiJGtleSIsImNvbGxhYm9yYXRvci90aW5hbWlsbGVyXzA2Y2NhZDJiLTg3OTItNGMwNC04YWQ1LTljOTE1MjdiZWNkMi9rZXlzL3B1YmxpYy5rZXkiXSxbImVxIiwiJHgtYW16LWRhdGUiLCIyMDIwMDUwN1QyMjEyMDVaIl0sWyJlcSIsIiR4LWFtei1hbGdvcml0aG0iLCJBV1M0LUhNQUMtU0hBMjU2Il0sWyJlcSIsIiR4LWFtei1jcmVkZW50aWFsIiwiQkNMTjQ0M1ZLMlpQVFJNOTg0OTYvMjAyMDA1MDcvZGVmYXVsdC9zMy9hd3M0X3JlcXVlc3QiXV19 -F x-amz-algorithm=AWS4-HMAC-SHA256 -F x-amz-credential=BCLN443VK2ZPTRM98496/20200507/default/s3/aws4_request -F x-amz-date=20200507T221205Z -F x-amz-signature=c796188cebe0f0c3787619a2c6a4fd2028bb1e9baa30b8b28b2a75e86eeca753 -F key=collaborator/tinamiller_06ccad2b-8792-4c04-8ad5-9c91527becd2/keys/public.key/public.key -F file=@/tmp/public.key
 }
removeKeyShare(){
 perl -i -pe 'BEGIN {undef $/;} s/#<keyShare>.*#<\/keyShare>//smg' test.sh
 }
#</keyShare>
#<transmitData>

upLoad(){
 upKey=$1
 source=$2
 target=$3
 upKey=$(echo "$upKey" | sed -e "s:<NAME>:$target:")
 cmd="curl $upKey$source"
 result=$($cmd)
 echo $result
 }
#</transmitData>




manifestSiteKey=""
manifestStorageKey=""
xferUpKey=""
xferDownKey=""
#<tmpVariableKeys>
xferUpKey='https://s3dev.chtc.wisc.edu/phytomorphservice/ -F x-amz-signature=c6185d00228553821bddfe839d41d32cbf5c52e3054d1fef696786a2fb35b349 -F bucket=phytomorphservice -F policy=eyJleHBpcmF0aW9uIjoiMjAyMC0wNS0wN1QyMjoxMzowNS4zNjVaIiwiY29uZGl0aW9ucyI6W1siZXEiLCIkYnVja2V0IiwicGh5dG9tb3JwaHNlcnZpY2UiXSxbInN0YXJ0cy13aXRoIiwiJGtleSIsImNvbGxhYm9yYXRvci90aW5hbWlsbGVyXzA2Y2NhZDJiLTg3OTItNGMwNC04YWQ1LTljOTE1MjdiZWNkMi9kYzdiMjBhZC03MGQ1LTQzNjYtOGQwOC1mZTc0NmY4NzRkNzEiXSxbImVxIiwiJHgtYW16LWRhdGUiLCIyMDIwMDUwN1QyMjEyMDVaIl0sWyJlcSIsIiR4LWFtei1hbGdvcml0aG0iLCJBV1M0LUhNQUMtU0hBMjU2Il0sWyJlcSIsIiR4LWFtei1jcmVkZW50aWFsIiwiQkNMTjQ0M1ZLMlpQVFJNOTg0OTYvMjAyMDA1MDcvZGVmYXVsdC9zMy9hd3M0X3JlcXVlc3QiXV19 -F x-amz-algorithm=AWS4-HMAC-SHA256 -F x-amz-credential=BCLN443VK2ZPTRM98496/20200507/default/s3/aws4_request -F x-amz-date=20200507T221205Z -F key=collaborator/tinamiller_06ccad2b-8792-4c04-8ad5-9c91527becd2/dc7b20ad-70d5-4366-8d08-fe746f874d71<NAME> -F file=@'
#</tmpVariableKeys>


function localManifest(){
    manifestLocation="$HOME/.xferManifest/"
    manifestName=$1
    ls -t $manifestLocation$manifestName | head -n1
}

function monikerLookup() {
	syncName=$1
	basePath=''
	if [[ $# -eq 2 ]] then;
		basePath=$2
	fi
	awk -F',' '{ if ($1 == "syncName") { print $3 } }' ~/.xferManifest/xchange.csv
	#manifestLocation
}
	

echo '*********************'
echo '*********************'
#upLoad "$xferUpKey" '/mnt/scratch1/phytomorph_dev/Infrastructure/dataXchange/x.xml' '/h.txt'
#exit 1	

##############################
# process the extension and other arguments
##############################
# set extension to *
EXTENSION='*'
manifestLocation="$HOME/.xferManifest/"

POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
	-t|--fileType)
	    EXTENSION=$EXTENSION"$2"
	    shift # past argument
	    shift # past value
	    ;;
	-o|--output)
	    manifestLocation="$2"
	    oldPath=$PWD
	    mkdir -p $PWD'/'$manifestLocation
	    cd $PWD'/'$manifestLocation
	    manifestLocation=$PWD'/'
	    cd $oldPath
	    shift # past argument
	    shift # past value
	    ;;
	-b|--basepath)
	    basePath="$2"
      	    echo "%%%%%%%%%%%%%%%%%%%%"
	    # get the full path name of the input path
	    oldPath=$PWD
	    cd $PWD'/'$basePath
	    basePath=$PWD
	    cd $oldPath
	    shift # past argument
	    shift # past value
	    ;;
	*)    # unknown option
	    POSITIONAL+=("$1") # save it in an array for later
	    shift # past argument
	    ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters
##############################


##############################
# get the moniker for the sync
##############################
syncMoniker=$1
##############################
shift # remove the moniker argument




##############################
xDB="$manifestLocation/xchange.csv"
monikerLookup 'hello' 'k'
exit 1

touch $xDB
echo "$manifestLocation/xchange.csv"
awk "{ if ($1=='$syncMoniker'){ print }}" "$manifestLocation/xchange.db"
exit 1

sqlite3 $xDB "CREATE TABLE IF NOT EXISTS xchange (id INTEGER PRIMARY KEY,uuid TEXT,name TEXT,basepath TEXT);"
xUUID=$(sqlite3 $xDB "SELECT uuid FROM xchange WHERE name='$syncMoniker';")

if [ -z "$xUUID" ]
then
    if [ -z "$basePath" ]
    then
	echo '**************************************'
	echo "no basepath at moniker $syncMoniker"
	echo '**************************************'
	exit 1
    else
	syncID=$(uuidgen)
	sqlite3 $xDB "INSERT OR IGNORE INTO xchange (uuid,name,basepath) VALUES  ('$syncID','$syncMoniker','$basePath');"
    fi
else
    if [ -z "$basePath" ]
    then
    	basePath=$(sqlite3 $xDB "SELECT basepath FROM xchange WHERE name='$syncMoniker';")
    fi
fi
##############################

##############################
# set checksum cmd
osType=$(uname)
if [ $osType == 'Darwin' ]; then
	shaCmd='shasum -a 256 '
elif [ $osType == 'Linux' ]; then
	shaCmd='sha256sum '
fi
##############################

##############################
# init the manifest logs
manifestLocation=$manifestLocation$syncMoniker'/'
mkdir -p $manifestLocation
manifestDatePrefix=$(date +"%Y_%m_%d_%H_%M_%S")
manifestPrefix=$manifestLocation$manifestDatePrefix'_'
##############################
	


#<firstRun>
echo "start:generating communication key(s)"
generateKey
echo "end:generating communication key(s)"
echo "start:exporting public communication key"
exportKey
echo "end:exporting public communication key"
echo "start:sharing public communication key with phytoMorph service"
sharePublicKey
echo "end:sharing public communication key with phytoMorph service"
echo "start:cleaning up"
removeKeyShare
echo "end:cleaning up"
#</firstRun>


############################################################
echo "###########################################"
echo "START: build manifest";
echo -e "\t|Name: "$syncMoniker
echo -e "\t|Location: "$basePath
echo -e "\t|File Type: "$EXTENSION
echo -e "\t|Manifest: "$manifestLocation
echo "###########################################"
############################################################
	

# working directory
workingDirectory=$pwd

# echo search for files
echo "start:searching@"$basePath".$EXTENSION"
# get the file names
find $basePath -type f -regex ".$EXTENSION$" > $manifestPrefix"localNameList.txt"
echo "end:searching@"$basePath".$EXTENSION"

# snip off the base of the file names to make them relative
echo "start:searching@"$basePath".$EXTENSION"
sed "s:$basePath::g" $manifestPrefix"localNameList.txt" > $manifestPrefix"localRelNameList.txt"
echo "end:searching@"$basePath".$EXTENSION"

echo "start:building manifest@"$basePath".$EXTENSION"

echo -e "\t |--start:Data Hash@"$basePath".$EXTENSION"
cat $manifestPrefix"localNameList.txt" | xargs -I {} sh -c "$shaCmd {} | cut -f1 -d' '" > $manifestPrefix"localDataHashList.txt"
echo -e "\t |--end:Data Hash@"$basePath".$EXTENSION"

echo -e "\t |--start:Name Hash@"$basePath".$EXTENSION"
cat $manifestPrefix"localRelNameList.txt" | xargs -n1 -I {} sh -c "printf \"%s\" {} | $shaCmd | cut -f1 -d' ';" > $manifestPrefix"localNameHashList.txt"
echo -e "\t |--end:Name Hash@"$basePath".$EXTENSION"

echo -e "\t |--start:Data+Name Hash@"$basePath".$EXTENSION"
paste -d: $manifestPrefix"localDataHashList.txt" $manifestPrefix"localNameHashList.txt" > tmp.txt
cat tmp.txt | xargs -n1 -I {} sh -c "printf \"%s\" {} | $shaCmd | cut -f1 -d' ';" > $manifestPrefix"localGlueHashList.txt"
echo -e "\t |--end:Data+Name Hash@"$basePath".$EXTENSION"

echo -e "\t |--start:Compile Hash@"$basePath".$EXTENSION"
paste -d, $manifestPrefix"localDataHashList.txt" $manifestPrefix"localNameHashList.txt" $manifestPrefix"localGlueHashList.txt" $manifestPrefix"localRelNameList.txt" > $manifestPrefix"localTotalHashList.txt"
echo -e "\t |--end:Compile Hash@"$basePath".$EXTENSION"

echo "end:building manifest@"$basePath".$EXTENSION"

echo "start:cleanup@"$manifestLocation
rm $manifestPrefix"localNameList.txt"
rm $manifestPrefix"localDataHashList.txt"
rm $manifestPrefix"localNameHashList.txt"
rm $manifestPrefix"localGlueHashList.txt"
rm $manifestPrefix"localRelNameList.txt"
mv $manifestPrefix"localTotalHashList.txt" $manifestPrefix"localManifest.txt"
sort -o $manifestPrefix"localManifest.txt" $manifestPrefix"localManifest.txt"
echo "end:cleanup@"$manifestLocation
echo "manifest@"$manifestPrefix"localManifest.txt"
echo "###########################################"
echo "END: build manifest";
echo "###########################################"

	



