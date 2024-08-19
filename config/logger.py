#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2018-2024 Institute of Computer Science of the Czech Academy of
# Sciences, Prague, Czech Republic. Authors: Martin Bures, Jaroslav Resler.
#
# This file is part of PALM-GeM.
#
# PALM-GeM is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# PALM-GeM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PALM-GeM. If not, see <https://www.gnu.org/licenses/>.

import sys
import logging

__all__ = ['extra_verbose', 'sql_extra_verbose', 'sql_debug', 'sql_verbose', 'verbose', 'debug',
           'warning', 'error', 'progress', 'change_log_level', 'logging_level', 'restore_log_level', 'logging_level_wr']

Logger = logging.getLogger('Logger')

EXTRA_VERBOSE = 1
def extra_verbose(msg, *args, **kwargs):
    source = {'source': 'EXTRA VERBOSE'}
    if logging.getLogger().isEnabledFor(EXTRA_VERBOSE):
        if args or kwargs:
            logging.log(EXTRA_VERBOSE, msg.format(*args, **kwargs), extra=source)
        else:
            logging.log(EXTRA_VERBOSE, msg, extra=source)

logging.addLevelName(EXTRA_VERBOSE, "EXTRA_VERBOSE")
logging.extra_verbose = extra_verbose

SQL_EXTRA_VERBOSE = 1
def sql_extra_verbose(con, *args, **kwargs):
    source = {'source': 'SQL EXTRA VERBOSE'}
    if logging.getLogger().isEnabledFor(SQL_EXTRA_VERBOSE):
        for n in con.notices:
            for nn in n.splitlines():
                logging.log(SQL_EXTRA_VERBOSE, nn.strip('NOTICE:').lstrip(), extra=source)
        con.notices.clear()

logging.addLevelName(SQL_EXTRA_VERBOSE, "SQL_EXTRA_VERBOSE")
logging.sql_extra_verbose = sql_extra_verbose

VERBOSE = 2
def verbose(msg, *args, **kwargs):
    source = {'source': 'VERBOSE'}
    if logging.getLogger().isEnabledFor(VERBOSE):
        if args or kwargs:
            logging.log(VERBOSE, msg.format(*args, **kwargs), extra=source)
        else:
            logging.log(VERBOSE, msg, extra=source)

logging.addLevelName(VERBOSE, "VERBOSE")
logging.verbose = verbose

SQL_VERBOSE = 2
def sql_verbose(con, *args, **kwargs):
    source = {'source': 'SQL VERBOSE'}
    if logging.getLogger().isEnabledFor(SQL_VERBOSE):
        for n in con.notices:
            for nn in n.splitlines():
                logging.log(SQL_VERBOSE, nn.strip('NOTICE:').lstrip(), extra=source)
        con.notices.clear()

logging.addLevelName(SQL_VERBOSE, "SQL_VERBOSE")
logging.sql_verbose = sql_verbose

DEBUG = 3
def debug(msg, *args, **kwargs):
    source = {'source': 'DEBUG'}
    if logging.getLogger().isEnabledFor(DEBUG):
        if args or kwargs:
            logging.log(DEBUG, msg.format(*args, **kwargs), extra=source)
        else:
            logging.log(DEBUG, msg, extra=source)

logging.addLevelName(DEBUG, "DEBUG")
logging.debug = debug

SQL_DEBUG = 3
def sql_debug(con):
    source = {'source': 'SQL DEBUG'}
    if logging.getLogger().isEnabledFor(SQL_DEBUG):
        for n in con.notices:
            for nn in n.splitlines():
                logging.log(SQL_DEBUG, nn.strip('NOTICE:').lstrip(), extra=source)
        con.notices.clear()

logging.addLevelName(SQL_DEBUG, "DEBUG")
logging.sql_debug = sql_debug

PROGRESS = 4
def progress(msg, *args, **kwargs):
    source = {'source': 'PROGRESS'}
    if logging.getLogger().isEnabledFor(PROGRESS):
        if args or kwargs:
            logging.log(PROGRESS, msg.format(*args, **kwargs), extra=source)
        else:
            logging.log(PROGRESS, msg, extra=source)

logging.addLevelName(PROGRESS, "PROGRESS")
logging.progress = progress


WARNING = 5
def warning(msg, *args, **kwargs):
    source = {'source': 'WARNING'}
    if logging.getLogger().isEnabledFor(WARNING):
        if args or kwargs:
            logging.log(WARNING, msg.format(*args, **kwargs), extra=source)
        else:
            logging.log(WARNING, msg, extra=source)

logging.addLevelName(WARNING, "WARNING")
logging.warning = warning

ERROR = 6
def error(msg, *args, **kwargs):
    source = {'source': 'ERROR'}
    if logging.getLogger().isEnabledFor(ERROR):
        if args or kwargs:
            logging.log(ERROR, msg.format(*args, **kwargs), extra=source)
        else:
            logging.log(ERROR, msg, extra=source)

logging.addLevelName(ERROR, "ERROR")
logging.error = error

Log_Format = "{asctime:20s} {source:15s} - {message}"

class LoggerWriter(object):
    def __init__(self, writer):
        self._writer = writer
        self._msg = ''

    def write(self, message):
        self._msg = self._msg + message
        while '\n' in self._msg:
            pos = self._msg.find('\n')
            self._writer(self._msg[:pos])
            self._msg = self._msg[pos+1:]

    def flush(self):
        if self._msg != '':
            self._writer(self._msg)
            self._msg = ''

# Log_Format = "%(levelname)-15s %(asctime)-5s %(source)-8s - %(message)s"
def logging_level(cfg):
    """ wrapper """
    logging_level_wr(cfg.logs.level, cfg.logs.path)

def logging_level_wr(level, path):
    logging.basicConfig(level=int(level),
                        format=Log_Format,
                        datefmt='%Y-%m-%d %H:%M:%S',
                        style='{',
                        handlers=[
                            logging.FileHandler(path,"w"),
                            logging.StreamHandler(sys.stdout),
                                 ])

    # sys.stdout = LoggerWriter(debug)
    sys.stderr = LoggerWriter(error)

    progress('Level of logging is {}', level)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)

def change_log_level(new_level):
    logging.getLogger().setLevel(new_level)
    progress('Level of logs was changed into level {}', new_level, logging._levelToName[new_level])

def restore_log_level(cfg):
    logging.getLogger().setLevel(cfg.logs.level)
    progress('Level of logs was changed back into level {}', cfg.logs.level, logging._levelToName[cfg.logs.level])
