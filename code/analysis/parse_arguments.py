#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 10:06:03 2022

@author: maxthomas_UO

parse arguments
"""

import ast
import datetime
import argparse

def make_args(parser):
    parser.add_argument('--suite', type=comma_split_str, required=False, default='ci501')
    parser.add_argument('--variable', type=comma_split_str, required=False, default='sos')
    parser.add_argument('--dates', type=comma_split_str, required=False, default='20150101,20190116')
    parser.add_argument('--region', type=comma_split_str, required=False, default='circumpolar')
    parser.add_argument('--lats', type=comma_split_float, required=False, default='-90,-60')
    parser.add_argument('--lons', type=comma_split_float, required=False, default='160,162')
    parser.add_argument('--depths', type=comma_split_float, required=False, default='0')
    parser.add_argument('--diff_abs', type=comma_split_str, required=False, default='absolute')
    parser.add_argument('--save_results', type=ast.literal_eval, required=False, default=True)
    parser.add_argument('--save_data', type=ast.literal_eval, required=False, default=True)
    parser.add_argument('--plot_sia_contour', type=ast.literal_eval, required=False, default=False)
    parser.add_argument('--plot_shelf', type=ast.literal_eval, required=False, default=True)
    parser.add_argument('--plot_0_wind', type=ast.literal_eval, required=False, default=False)
    parser.add_argument('--land_mask', type=ast.literal_eval, required=False, default=False)
    parser.add_argument('--base_suite',type=str, required=False, default='cmip')
    parser.add_argument('--analysis',type=comma_split_str, required=False, default='timeseries')
    parser.add_argument('--experiment_id', type=str, required=False, default='ssp585')
    parser.add_argument('--season', type=str, required=False, default='year')
    parser.add_argument('--abs_lims', type=comma_split_float, required=False, default=-999)
    parser.add_argument('--diff_lims', type=comma_split_float, required=False, default=-999)
    parser.add_argument('--computer', type=str, required=False, default='maui')
    parser.add_argument('--run_type', type=str, required=False, default='main')
    args = vars(parser.parse_args())
    args['script'] = parser.prog
    args['run_time'] = datetime.datetime.now()
    return args

def comma_split_str(arg_in):
    return arg_in.split(',')

def comma_split_float(arg_in):
    return list(map(float, arg_in.split(',')))

def inputs2arguments():
    return make_args(argparse.ArgumentParser())



