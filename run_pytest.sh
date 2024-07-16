#!/bin/bash

/home/user/.local/bin/pytest

EXIT_CODE=$?
echo "pytest exited with code: $EXIT_CODE"
exit $EXIT_CODE
