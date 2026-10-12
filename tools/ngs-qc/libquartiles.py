#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on 27/03/15

@author: Matthias Blum
"""

import json
import os
import sqlite3
import sys


def init_db(dbjson):
    with open(dbjson) as fp:
        versions = json.load(fp)

    dbfile = os.path.join(os.path.dirname(__file__), 'db', 'quartiles.db')

    if os.path.isfile(dbfile):
        os.unlink(dbfile)


    con = sqlite3.connect(dbfile)
    cur = con.cursor()
    cur.execute("CREATE TABLE quartile ("
                "disp REAL NOT NULL, "
                "quartile INTEGER NOT NULL, "
                "value READ NOT NULL, "
                "version TEXT NOT NULL, "
                "active INTEGER NOT NULL)")

    i = 1
    for v, version in sorted(versions.items()):
        for disp, quartiles in version.items():
            for q, value in quartiles.items():
                cur.execute("INSERT INTO quartile (disp, quartile, value, version, active) "
                            "VALUES (?, ?, ?, ?, ?)", (disp, q, value, v, i == len(versions)))
        i += 1

    con.commit()
    cur.close()
    con.close()


def get():
    """
    Get the quartiles score for each dispersion of a given background subtraction mode
    :param bg:  "global" for Poisson distribution (only one currently supported)
    :return:
    """
    dbfile = os.path.join(os.path.dirname(__file__), 'db', 'quartiles.db')
    con = sqlite3.connect(dbfile)
    cur = con.cursor()
    cur.execute('SELECT disp, quartile, value, version FROM quartile WHERE active=1')
    res = cur.fetchall()
    cur.close()
    con.close()

    quartiles = {}
    version = None
    for row in res:
        version = row[3]
        if row[0] not in quartiles:
            quartiles[row[0]] = {row[1]: row[2]}
        else:
            quartiles[row[0]][row[1]] = row[2]

    return quartiles, version


if __name__ == '__main__':
    init_db(sys.argv[1])