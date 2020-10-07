#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 13:11:02 2016

@author: ioannis vogiatzis oikonomidis
         ioannis.vogiatzis@psi.ch
"""

import re
import xml.etree.ElementTree as ET
import sys
import os
import argparse

class logFileParser:

    '''
    class that imports an xml or log file and converts it into a dictionary
    '''

    def __init__(self):
        self.logDict = {}

    def xmlParser(self, filename):
        self.logDict = {}

        '''
        imports and xml file as a tree and then converts into a python dictionary
        by itetarating along all nodes of the tree and using their name as a dictionary key'''
        tree = ET.parse(filename)

        for element in tree.getiterator():
            if len(element.items()) or len(element.keys()):
                key = [el[1] for el in element.items()]
                key = ' '.join(key)

                #check if you have a number or a string
                value = element.text
                try:
                    value = float(value)
                except ValueError:
                    value
                #check if number is integer
                if type(value) == float:
                    value = str(value)
                    if len(value.split('.')[-1]) == 1 and value.split('.')[-1] == '0':
                        value = int(float(value))
                    else:
                        value = float(value)

                self.logDict.update({key.title().replace(' ', ''): value})
            else:
                try:
                    self.logDict.update({element.tag: float(element.text)})
                except:
                    self.logDict.update({element.tag: element.text})

            root = tree.getroot()
            children = root.getchildren()
            self.logDict['name'] = children[0].text

    def logParser(self, filename):
        self.logDict = {}

        prefix = ""
        with open(filename, "r") as cache:
            # read file into a list of lines
            lines = cache.readlines()
            # loop through lines
            for line in lines:
                # skip lines starting with "--".
                if line.startswith("--"):
                    #create a prefix for each dictionary entry from dashed lines
                    line = re.sub("-", "\t", line).strip().split('\t')
                    prefix = line[0]
                elif not ':' in line:
                    continue
                else:
                    # '\s*' matches any whitespace zero or more times
                    # ':\s' matches the colon and the trailing space.
                    # We thus split the line at 'space(s) colon space(s)' or
                    # 'colon space(s)'
                    # Strip the trailing return (\n)
                    # Split into list using '\t' as the split pattern
                    # Test RegEx with https://regex101.com/
                    line = re.sub("\s*:\s*", "\t", line).strip().split('\t')
                    key = prefix + " " + line[0]

                    # add prefix to the entry
                    # use first item in list for the key, join remaining list items
                    # with ", " for the value.
                    #logDict[key] = ", ".join(line[1:])
                    value = line[1:][0]
                    try:
                        value = float(value)
                    except ValueError:
                        value
                        #check if number is integer
                    if type(value) == float:
                        value = str(value)
                        if len(value.split('.')[-1]) == 1 and value.split('.')[-1] == '0':
                            value = int(float(value))
                        else:
                            value = float(value)
                    self.logDict[key] = value

#    def logParser(self, filename):
#        self.logDict = {}
#        '''
#        Fucnction written by David Haberthür and taken from
#        https://gist.github.com/habi/7a43054c9cf3f1357755#file-samplelogfile-log-L
#        '''
#        with open(filename, "r") as cache:
#            # read file into a list of lines
#            lines = cache.readlines()
#            # loop through lines
#            for line in lines:
#                # skip lines starting with "--".
#                if not line.startswith("--"):
#                    # '\s*' matches any whitespace zero or more times
#                    # ':\s' matches the colon and the trailing space.
#                    # We thus split the line at 'space(s) colon space(s)' or
#                    # 'colon space(s)'
#                    # Strip the trailing return (\n)
#                    # Split into list using '\t' as the split pattern
#                    # Test RegEx with https://regex101.com/
#                    line = re.sub("\s*:\s*", "\t", line).strip().split('\t')
#                    # use first item in list for the key, join remaining list items
#                    # with ", " for the value.
#                    self.logDict[line[0]] = ", ".join(line[1:])

    def searchValueOfKey(self, dictionary, lookup, verbose=False):
#        '''
#        Fucnction written by David Haberthür and taken from
#        https://gist.github.com/habi/7a43054c9cf3f1357755#file-samplelogfile-log-L
#        '''
#        if verbose:
#            print 'Looking for "%s"' % lookup
#        for key, value in dictionary.items():
#            if str(lookup.lower()) in str(key) or str(lookup.upper()) in str(key):
#                if verbose:
#                    print 'We found the "%s" to be "%s"' % (key, value)
#                return value


            #old function was case sentisive

        if verbose:
            print('Looking for "%s"' % lookup)
        for key, value in dictionary.items():
            if re.search(lookup, key.replace(' ', ''), re.IGNORECASE):
                if verbose:
                    print('We found the "%s" to be "%s"' % (key, value))
                return value

    def parser_input(self, filename):
        if filename.endswith('xml'):
            self.xmlParser(filename)
        elif filename.endswith('log'):
            self.logParser(filename)
        else:
            sys.exit('Logfile type given not suported. Exiting...')

def getArgs():
    argsparser = argparse.ArgumentParser(description='Import a log file and convert it into a python dictionary', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    argsparser.add_argument('-Di', '--pathin', dest='pathin', default='./', type=str, help='Specify input file')
    argsparser.add_argument('-v', '--verbose', help="increase output verbosity", action="store_true")

    args = argsparser.parse_args()

    ##  Exit of the program in case the compulsory arguments,
    ##  are not specified
    if args.pathin is None:
        argsparser.print_help()
        sys.exit('ERROR: Invalid log file!')

    return args

def main():
    args = getArgs()

    if not os.path.isfile(args.pathin):
        print("directory doesn't exist. exiting..")
        sys.exit()

    parser = logFileParser()

    parser.parser_input(args.pathin)

    if args.verbose == True:
        print(parser.logDict)

if __name__ == '__main__':

    main()
	