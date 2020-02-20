
OLD=$(cat ./pilib/commit_info.txt)
NEW="const std::string VERSION = \"$(git describe --dirty --always --tags)\";"

#echo $OLD
#echo $NEW

# Only write the file if the version is different.
# This is to ensure that no re-build is triggered unnecessarily.
if [ "$OLD" != "$NEW" ]
then
    echo $NEW > ./pilib/commit_info.txt
fi

