#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#Marcelo D. Poleto
#MAR 2023


import sys
from configparser import ConfigParser
from configparser import NoOptionError
import numpy as np


##############################################################
# Parsing a configuration file and returning a Config object
#
# - configuration parameters are returned as object attributes
# - exception errors are thrown for wrong values or wrong
#   section headers are provided
#
##############################################################

class Configuration:
    def __init__(self, configinput):
        self.config = ConfigParser(allow_no_value=True, inline_comment_prefixes="#")
        self.configinput = configinput
        self.parse_config()


    def parse_config(self):
        if self.configinput is None:
            raise ValueError(">>> Configuration file is a necessary input. Use -template to generate one.\n>>> Exiting...\n")
        else:
            self.config.read(self.configinput)

        allowed_sections = ["Environment Selection", "Probe Selection", "Solvent", "Time", "Box Info"]
        for section in self.config.sections():
            if section not in allowed_sections:
                raise ValueError(""">>> ERROR: Unknown section ('""" + section + """') in configuration file!\n""")

        self.parse_environment_selection()
        self.parse_probe_selection()
        self.parse_solvent()
        self.parse_time()
        self.parse_box_info()


    def parse_environment_selection(self):
        try:
            self.sele_elecfield = self.config['Environment Selection']["sele_environment"].strip('"')
        except:
            raise ValueError("\n>>> ERROR: sele_environment must be defined with a valid MDanalysis selection!\n")

    def parse_probe_selection(self):
        try:
            self.mode = self.config['Probe Selection']["mode"].strip('"').lower()
            if self.mode == "atom":
                try:
                    self.selatom = self.config['Probe Selection']["selatom"].strip('"')
                except:
                    raise ValueError("""\n>>> ERROR: in "ATOM" mode, selatom must be defined!\n""")
            elif self.mode == "bond":
                try:
                    self.selbond1 = self.config['Probe Selection']["selbond1"].strip('"')
                    self.selbond2 = self.config['Probe Selection']["selbond2"].strip('"')
                except:
                    raise ValueError(""">>> ERROR: in "BOND" mode, both selbond1 and selbond2 must be defined!\n""")
            elif self.mode == "coordinate":
                try:
                    tmp_probecoordinate = self.config['Probe Selection']["probecoordinate"].strip('[]').split(",")
                    self.probecoordinate = np.array([float(item) for item in tmp_probecoordinate])
                except:
                    raise ValueError(""">>> ERROR: in "COORDINATE" mode, a list of coordinates [X,Y,Z] must be provided!\n""")
            elif self.mode == "list":
                try:
                    self.file_of_coordinates = self.config['Probe Selection']["file_of_coordinates"].strip('"')
                except:
                    raise ValueError(""">>> ERROR: in "LIST" mode, "file_of_coordinates" must be defined!\n""")
            else:
                raise ValueError(""">>> ERROR: "mode" must be defined as "ATOM", "BOND", "COORDINATE" or "LIST"!\n""")
        except:
            raise ValueError(""">>> ERROR: "mode" must be defined as "ATOM", "BOND", "COORDINATE" or "LIST"!\n""")

        if self.mode == "coordinate" or self.mode == "list":
            try:
                self.remove_self = self.config['Probe Selection'].getboolean("remove_self")
                if self.remove_self:
                    try:
                        self.remove_cutoff = float(self.config['Probe Selection']["remove_cutoff"])
                    except KeyError:
                        print("""\n>>> WARNING: "remove_cutoff" was not provided! Using the default value (1 Angstrom)!...\n""")
                        self.remove_cutoff = 1
                else:
                    self.remove_cutoff = 0
            except:
                raise ValueError(""">>> ERROR: "remove_self" must be a boolean!\n)""")
                self.remove_self = False
                self.remove_cutoff = 0
        else:
            self.remove_self = False
            self.remove_cutoff = 0


    def parse_solvent(self):
        try:
            self.include_solvent = str(self.config['Solvent']["include_solvent"].strip('"')).lower()
            if self.include_solvent == "true":
                self.include_solvent = True
                try:
                    self.solvent_selection = str(self.config['Solvent']["solvent_selection"].strip('"'))
                except:
                    raise ValueError(""">>> ERROR: solvent_selection must be a valid MDanalysis selection!\n""")
                try:
                    self.solvent_cutoff = float(self.config['Solvent']["solvent_cutoff"])
                except:
                    raise ValueError(""">>> ERROR: solvent_cutoff must be a number!\n""")
            else:
                self.include_solvent = False
        except:
            self.include_solvent = False

    def parse_time(self):
        try:
            self.dt = int(self.config['Time']['dt'])
        except:
            self.dt = 1


    def parse_box_info(self):
        try:
            if 'Box Info' in self.config.sections() and 'redefine_box' in self.config['Box Info']:
                redefine_box = str(self.config['Box Info']['redefine_box'].strip('"')).lower()
                self.redefine_box = True if redefine_box == 'true' else False

                if self.redefine_box:
                    boxdimensions_str = str(self.config['Box Info']['boxdimensions'].strip('][')).lower()
                    self.boxdimensions = [float(item) for item in boxdimensions_str.split(',')]

                    if len(self.boxdimensions) != 6:
                        raise ValueError('boxdimensions should contain [a,b,c, alpha, beta, gamma]!')
                else:
                    self.redefine_box = False
                    self.boxdimensions = None
            else:
                self.redefine_box = False
                self.boxdimensions = None
        except:
            raise ValueError('To redefine the box, boxdimensions must be provided!')
