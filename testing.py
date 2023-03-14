import difflib
import os
from glob import glob
import re

HEADER = '\033[95m'
OKBLUE = '\033[94m'
OKCYAN = '\033[96m'
OKGREEN = '\033[92m'
OKRED = '\033[33m'
WARNING = '\033[93m'
FAIL = '\033[91m'
ENDC = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'

FILE_TO_TEST = "tsp"

PATH_TO_TESTS = "pub-instances/*"


def PrintFail(arg, name):
    if arg is False:
        print(f"{OKRED}{BOLD}\t-->\t{name.upper()} FAILED{ENDC}\n\n")
    return True

def PrintPass(name):
    print(f'{OKGREEN}{BOLD}\t-->\t{name.upper()} PASSED{ENDC}\n\n')


print("\n\t -----------------",
      f"\n\t{BOLD}{UNDERLINE}TESTING {FILE_TO_TEST}{ENDC}\n",
      "\t-----------------\n\n")

test_statistics = {"total": 0, "failed": 0, "passed": 0}

test_files = sorted(glob(PATH_TO_TESTS))
for test in test_files:
    test_statistics["total"] += 1

    _, name = os.path.split(test)

    print(f'  {OKBLUE}{name.upper()} --| {ENDC}')

    r1 = re.findall(r"\d+", test)

    os.system(f'./{FILE_TO_TEST} {test} {r1[-1]} > our_output/{test}')

    output = open(f"output/{test}").readlines()
    tmp = open(f'our_output/{test}').readlines()

    thereIsDiff = False
    for line in difflib.unified_diff(output, tmp):
        thereIsDiff = PrintFail(thereIsDiff, name)
        print(line)
    if not thereIsDiff:
        PrintPass(name)
        test_statistics["passed"] += 1
    else:
        test_statistics["failed"] += 1
    
if test_statistics["failed"] > 0:
    print(
        f'{BOLD}{OKRED}\t FAILED -> TOTAL: {test_statistics["total"]} | PASSED: {test_statistics["passed"]} | FAILED: {test_statistics["failed"]}{ENDC}')
else:
    print(
        f'{BOLD}{OKGREEN}\t PASSED -> TOTAL: {test_statistics["total"]} | PASSED: {test_statistics["passed"]} | FAILED: {test_statistics["failed"]}{ENDC}')

print()
