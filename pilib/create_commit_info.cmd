@echo off
SetLocal EnableDelayedExpansion

IF NOT EXIST commit_info.txt (
	set OLD=
) ELSE (
	set /p OLD=<commit_info.txt
)


set NEW=unknown

where git >nul 2>nul
IF !ERRORLEVEL! EQU 0 (
	git rev-parse --is-inside-work-tree >nul 2>nul
    IF !ERRORLEVEL! EQU 0 (
		git describe --dirty --always --tags > commit_info_temp.txt
		set /p NEW=<commit_info_temp.txt
	)
)

echo std::string VERSION = ^"!NEW!^"; > commit_info_temp.txt
set /p NEW=<commit_info_temp.txt

rem echo %OLD%
rem echo %NEW%

IF NOT !OLD!==!NEW! (
	echo Version string changed to !NEW!
	echo !NEW!>commit_info.txt
)
