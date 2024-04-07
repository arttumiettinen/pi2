
# Check if git is installed and if the source is in a git repo.
if command -v git &> /dev/null
then

	if git rev-parse --is-inside-work-tree &> /dev/null
	then
		OLD=$(cat ./pilib/commit_info.txt 2>/dev/null)
		NEW="const std::string VERSION = \"$(git describe --dirty --always --tags)\";"

		#echo $OLD
		#echo $NEW

		# Only write the file if the version is different.
		# This is to ensure that no re-build is triggered unnecessarily.
		if [ "$OLD" != "$NEW" ]
		then
			echo $NEW > ./pilib/commit_info.txt
		fi
		
		exit
	fi
fi

# Git is not available or the source is not in a git repo.
# Output dummy version file if no version file has been provided.
if [ ! -f ./pilib/commit_info.txt ]; then
	echo "const std::string VERSION = \"unknown\";" > ./pilib/commit_info.txt
fi
