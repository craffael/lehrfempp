#!/usr/bin/env python3
# coding=utf-8
# Author: Huang Liaowang
# Date: December 2019
# This is part of the LehrFEM++ code suite

import os
import re
import sys
import warnings

"""
This script is to replaces all the @lref{<label>} statements in a C++ files with a label-specific 
string and the corresponding number from .aux file. Here we assume the line corresponding to this label in .aux file 
looks like

\newlabel{<label>@cref}{{[<Label>][x][xxx]<number>}{[x][xx][x]xxx}}

We would like to replace @lref{<label>} with <Label> <number>, which means that label is
extracted from the first square brackets and the number comes from the content after two
brackets of <Label>.  
"""

""" 
The ultimate goal is to have references to LaTeX documents in documentation generated
by Doxygen, and this script is to do the preprocessing before the Doxygen parse the source
file.

To use this filter, specifiy 
    INPUT_FILTER  = "python" "filter.py"
and put this script and "NCSE19refs.aux" together with doxygen configuration file.
"""

label_dict = {
    "equation": "Equation",
    "par": "Paragraph",
    "chapter": "Chapter",
    "sec": "Section",
    "subsection": "Subsection",
    "figure": "Figure",
    "code": "Code",
    "remark": "Remark"
}


def label_number(ref, label_dict):
    """

    :param ref: a line in .aux file, e.g. '\newlabel{eq:vec@cref}{{[equation][5][1,2,3]1.2.3.5}{[1][44][]44}}'
    :param label_dict: a dictionary, e.g. convert 'equation' to 'Equation'
    :return: label+number, e.g. Equation 1.2.3.4
    """
    ref_extract = re.split('[{}]', ref)[4]  # [equation][5][1,2,3]1.2.3.5
    ref_split = re.split('[\[\]]', ref_extract)  # ['', 'equation', '', '5', '', '1,2,3', '1.2.3.5']
    label = ref_split[1]
    number = ref_split[-1]
    if label in label_dict:
        label = label_dict[label]
    return label + ' ' + number


def aux_dict(aux_file):
    """
    parse the .aux file and store in a map.
    The line in .aux file we are interested is assumed to look like:
    \newlabel{<label>@cref}{{[<Label>][x][xxx]<number>}{[x][xx][x]xxx}}
    And we search for "@cref" to find the line.

    :param aux_file: the .aux file
    :return: a dictionary maps all the reference to label + number
    """
    auxdict = {}
    aux = open(aux_file, "r")
    for line in aux.readlines():
        ref = re.findall('{(.*?)@cref}', line)
        if len(ref) > 0:  # find the line which contains {xxx@cref}
            auxdict[ref[0]] = label_number(line, label_dict)
    return auxdict


def find_label_number(cref):
    """
    The line in .aux file we are searching for is assumed to look like:
    \newlabel{<label>@cref}{{[<Label>][x][xxx]<number>}{[x][xx][x]xxx}}
    And we search for "<lable>@cref" to find the line.

    :param cref: re match object, e.g. cref.group(0)='@lref{eq:vec}'
    :return:
    """
    ref = cref.group(0)
    if re.search(r'\n', ref, flags=re.DOTALL) is not None:
        warnings.warn("{} crosses the lines".format(ref), Warning)
        return ref
    clabel = re.findall('{(.*)}', ref)[0]  # get the content in {}, e.g. 'eq:vec'
    if clabel in auxdict:
        return auxdict[clabel]
    else:
        # raise a warning
        warning_str = "Cannot find a label for {}".format(ref)
        warnings.warn(warning_str, Warning)
        return ref
    # clabel = re.findall('{(.*)}', cref.group(0))[0]  # get the content in {}, e.g. 'eq:vec'
    # aux = open(aux_file, "r")
    # ref = None
    # for line in aux.readlines():
    #     if re.search(clabel + '@cref', line) is not None:  # find the line which has eq:vec@cref
    #         ref = line
    #         break
    # if ref is None:
    #     raise ValueError("{} not found in .aux file".format(clabel))
    # return label_number(ref, label_dict)


def replace_lref(file):
    """
    Replace the @lref{...} with specific label and number in a file, and print the result (stdout)

    :param file: the file needs to parse
    """

    # print("{} will be processed".format(file))

    # read file as string
    file_handle = open(file, 'r')
    file_string = file_handle.read()
    file_handle.close()

    # do the replacement
    # @lref.*?{.*?} can find @lref{} even they cross the lines
    pattern = re.compile(r'@lref{.*?}|@lref *\n// *?{.*?}|@lref *\n *?{.*?}', re.DOTALL)
    file_string = re.sub(pattern, find_label_number, file_string)
    # <.*?> will match only '<a>'. <.*> is greedy, it may match '<a> b <c>'

    print(file_string)  # to stdout


if __name__ == "__main__":
    # takes a file as input and print the content after parsing.
    # the output is on stdout which can be understood natively by doxygen.

    if len(sys.argv) == 2:
        aux_file = "NPDEFLrefs.aux"  # default .aux file
        file = sys.argv[1]
    else:
        aux_file = sys.argv[1]  # input argument
        file = sys.argv[2]
    auxdict = aux_dict(aux_file)
    replace_lref(file)