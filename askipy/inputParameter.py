# ----------------------------------------------------------------------------
#   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
#
#   This file is part of ASKI version 1.2.
#
#   ASKI version 1.2 is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.
#
#   ASKI version 1.2 is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
# ----------------------------------------------------------------------------
class inputParameter:
    """
    """
    def __init__(self, filename):
        self._key_val = dict([
                (line.split('=')[0].strip(),  line.split('=')[1].split('#')[0].strip())
                for line in open(filename, 'r')
                if line.strip() != ''            # ignore empty lines or lines containing spaces only
                if line.strip()[0:1] != '#'      # ignore comment lines (first non space character is '#')
                if '=' in line.split('#')[0]])   # ignore lines without '=' in front of '#' (works if no '#' in line)

    def hasKey(self, key):
        return key in self._key_val

    def keysNotPresent(self, keys):
        lnp = []
        for key in keys:
            if key not in self._key_val:
                lnp.append(key)
        return lnp

    def prnt(self):
        for item in self._key_val.items():
            print("'"+"' = '".join(item)+"'")

    def log(self, logfile):
        with open(logfile, "a") as fl:
            for item in self._key_val.items():
                fl.write("'"+"' = '".join(item)+"'"+"\n")

    def getval(self, key):
        if key in self._key_val:
            return self._key_val[key]
        else:
            return None

    def ival(self, key):
        if key in self._key_val:
            try:
                return int(self._key_val[key])
            except TypeError:
                return None
        else:
            return None
        
    def ilist(self, key, n=None):
        if key in self._key_val:
            try:
                li = [int(i) for i in self._key_val[key].split()]
                if not n:
                    return li
                if n > 0:
                    return li[:min(n, len(li))]
                else:
                    return None
            except TypeError:
                return None
        else:
            return None
        
    def fval(self, key):
        if key in self._key_val:
            try:
                return float(self._key_val[key])
            except TypeError:
                return None
        else:
            return None
        
    def flist(self, key, n=None):
        if key in self._key_val:
            try:
                lf = [float(i) for i in self._key_val[key].split()]
                if not n:
                    return lf
                if n > 0:
                    return lf[:min(n, len(lf))]
                else:
                    return None
            except TypeError:
                return None
        else:
            return None
        
    def sval(self, key):
        if key in self._key_val:
            try:
                return self._key_val[key]
            except TypeError:
                return None
        else:
            return None

    def slist(self, key, n=None):
        if key in self._key_val:
            try:
                ls = self._key_val[key].split()
                if not n:
                    return ls
                if n > 0:
                    return ls[:min(n, len(ls))]
                else:
                    return None
            except TypeError:
                return None
        else:
            return None
        
    def lval(self, key):
        if key in self._key_val:
            if self._key_val[key].lower() in ['.true.', 'true', 't', '1', '1.']:
                return True
            elif self._key_val[key].lower() in ['.false.', 'false', 'f', '0', '0.']:
                return False
            else:
                return None
        else:
            return None
