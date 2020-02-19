@echo off
SetLocal EnableDelayedExpansion

IF NOT EXIST commit_info.txt (
	set /p OLD=
) ELSE (
	set /p OLD=<commit_info.txt
)


git describe --dirty --always --tags > commit_info_temp.txt
set /p NEW=<commit_info_temp.txt
echo std::string VERSION = ^"%NEW%^"; > commit_info_temp.txt
set /p NEW=<commit_info_temp.txt

REM echo %OLD%
REM echo %NEW%

IF NOT !OLD!==!NEW! (
	echo Version string changed
	echo %NEW%>commit_info.txt
)

