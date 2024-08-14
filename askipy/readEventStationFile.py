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
class eventList:
    """
    This class reads event list files in analogy with module readEventStationFile.f90. Every specific event
    list file type supported in fortran module readEventStationFile should also be supported here.
    Supported list_type values are 'standard' and None
    """

    def __init__(self, filename, list_type=None):
        self.events = dict()
        self.nev = 0
        self.csys = ''
        if list_type is None or str(list_type).lower() == 'standard':
            self.create_standard(filename)
        else:
            raise Exception("event list type '"+str(list_type)+"' is not supported")

# ----------------------------------------------------------------------------

    @classmethod
    def standard(cls, filename):
        return cls(filename, list_type='standard')

# ----------------------------------------------------------------------------

    def create_standard(self, filename):
        """
        This routine reads a standard event list file
        A standard event list evlist has the attribute events.
        evlist.events is a dictionary which can be accessed by
        eventID or event index (offset 0 index). The values
        of that dictionary themselves are dictionaries containing
        the event information, having keys 'evid', 'otime', 'slat',
        'slon', 'sdepth', 'mag','styp', 'force' (if styp==0) and 'momten' (if styp==1)
        """
        # read lines of event list file
        try:
            lines = open(filename, 'r').readlines()
        except OSError:
            raise Exception("could not read in lines of event list file '"+filename)

        # initalize empty lists for dictionary self.events
        events_keys = []
        events_vals = []

        # check coordinate system of events
        csys = lines[0].strip().split()[0]
        if csys not in ['C', 'S']:
            raise Exception("first entry '"+csys+"' of first line of event list file '" + filename +
                            "' which defines the coordinate system is neither 'C' nor 'S'")

        # process valid lines from file
        for line in lines[1:]:
            vals = line.strip().split()

            if len(vals) < 8:
                continue

            # check type of source (force, moment tensor)
            try:
                styp = int(vals[7])
            except TypeError:
                continue
            if styp == 0:
                if len(vals) >= 11:
                    keys = ['evid', 'otime', 'slat', 'slon', 'sdepth', 'mag', 'tred', 'styp', 'force']
                    events_keys.append(vals[0])
                    events_vals.append(dict(zip(keys, vals[:8] + vals[8:11])))
                else:
                    continue
            elif styp == 1:
                if len(vals) >= 14:
                    keys = ['evid', 'otime', 'slat', 'slon', 'sdepth', 'mag', 'tred', 'styp', 'momten']
                    events_keys.append(vals[0])
                    events_vals.append(dict(zip(keys, vals[:8] + vals[8:14])))
                else:
                    continue
            else:
                keys = ['evid', 'otime', 'slat', 'slon', 'sdepth', 'mag', 'tred', 'styp']
                events_vals.append(dict(zip(keys, vals[:8])))

        # if there were any valid lines, finally define self.events
        nev = len(events_keys)
        if nev != 0:
            # self.events is a dictionary which can be accessed by eventID or event index (offset 0 index)
            # the values of dictionary self.events themselves are dictionaries containing the event information
            self.events = dict(zip(events_keys+list(range(nev)), 2*events_vals))
            self.nev = nev
            self.csys = csys


class stationList:
    """
    This class reads station list files in analogy with module readEventStationFile.f90.
    Every specific station list file type supported in fortran module readEventStationFile
    should also be supported here. Supported list_type values are 'standard' and None
    """

    def __init__(self, filename, list_type=None):
        self.stations = dict()
        self.nstat = 0
        self.csys = ''
        if list_type is None or str(list_type).lower() == 'standard':
            self.create_standard(filename)
        else:
            raise Exception("station list type '" + str(list_type) + "' is not supported")

    @classmethod
    def standard(cls, filename):
        return cls(filename, list_type='standard')

    def create_standard(self, filename):
        """
        This routine reads a standard station list file. A standard station list stlist has the attribute stations.
        stlist.stations is a dictionary which can be accessed by station name or station index (offset 0 index).
        The values of that dictionary themselves are dictionaries containing the station information,
        having keys 'staname', 'netcode', 'lat', 'lon' and 'alt'
        """
        try:
            lines = open(filename, 'r').readlines()
        except OSError:
            raise Exception("could not read in lines of station list file '" + filename)

        # initalize empty lists for dictionary self.stations
        stations_keys = []
        stations_vals = []

        # define keys for each station entry
        keys = ['staname', 'netcode', 'lat', 'lon', 'alt']

        # check coordinate system of events
        csys = lines[0].strip().split()[0]
        if csys not in ['C', 'S']:
            raise Exception("first entry '" + csys + "' of first line of station list file '" + filename +
                            "' which defines the coordinate system is neither 'C' nor 'S'")

        # process valid lines from file
        # choose netcode.staname as key for station to distinguish identical station names in different networks
        for line in lines[1:]:
            vals = line.strip().split()
            if len(vals) < 5:
                continue
            stations_keys.append(vals[1]+'.'+vals[0])              # set key for station
            stations_vals.append(dict(zip(keys, vals[:5])))        # set values for stations subdict

        # if there were any valid lines, finally define self.stations
        nstat = len(stations_keys)
        if nstat != 0:
            # self.stations is a dictionary which can be accessed by station name or station index (offset 0 index)
            # the values of dictionary self.stations themselves are dictionaries containing the station information
            self.stations = dict(zip(stations_keys+list(range(nstat)), 2*stations_vals))
            self.nstat = nstat
            self.csys = csys
