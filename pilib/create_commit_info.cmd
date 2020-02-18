
git describe --dirty --always --tags > commit_info.txt
set /p var=<commit_info.txt
echo std::string VERSION = ^"%var%^"; > commit_info.txt


