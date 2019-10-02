################################################################################
#
#   BackstopHistory - Class used for containing and assembling  a
#                     backstop history
#
#     Author: Gregg Germain
#
#   Original: BackstopHistory.py
#
#
################################################################################
from __future__ import print_function
import copy
import glob
import numpy as np
import os
import logging
from pathlib import Path
from importlib import reload

from backstop_history import LTCTI_RTS

import kadi.commands
from Chandra.Time import DateTime


# -------------------------------------------------------------------------------
#
#  globfile
#
# -------------------------------------------------------------------------------
def globfile(pathglob):
    """Return the one file name matching ``pathglob``.  Zero or multiple
    matches raises an IOError exception."""

    files = glob.glob(pathglob)
    if len(files) == 0:
        raise IOError('No files matching %s' % pathglob)
    elif len(files) > 1:
        raise IOError('Multiple files matching %s' % pathglob)
    else:
        return files[0]

def config_logger(verbose):
    """
    Set up console logger.

    Parameters
    ----------
    verbose : integer
        Indicate how verbose we want the logger to be.
        (0=quiet, 1=normal, 2=debug)
    """

    # Disable auto-configuration of root logger by adding a null handler.
    # This prevents other modules (e.g. Chandra.cmd_states) from generating
    # a streamhandler by just calling logging.info(..).
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass
    rootlogger = logging.getLogger()
    rootlogger.addHandler(NullHandler())
    logger = logging.getLogger("backstop_history")
    logger.setLevel(logging.DEBUG)

    # Set numerical values for the different log levels
    loglevel = {0: logging.CRITICAL,
                1: logging.INFO,
                2: logging.DEBUG}.get(verbose, logging.INFO)

    formatter = logging.Formatter('%(message)s')

    console = logging.StreamHandler()
    console.setFormatter(formatter)
    console.setLevel(loglevel)
    logger.addHandler(console)

    return logger


class BackstopHistory(object):

    def __init__(self, cont_file_name='ACIS-Continuity.txt',
                 NLET_tracking_file_path='/data/acis/LoadReviews/NonLoadTrackedEvents.txt',
                 logger=None, verbose=1):
        if logger is None:
            logger = config_logger(verbose)
        self.logger = logger

        # Full path to a test output file used for debugging
        self.debug_outfile_path = '/home/gregg/Desktop/Debug_cmds.dat'

        # Full list of backstop commands whihc have been assembled
        self.master_list = []

        self.rev_to_take = []

        # The time of the Review Load first command in DOY and seconds format
        self.Date_OFC = None
        self.Time_OFC = None

        # The time of the safemode/shutdown in DOY and seconds format
        self.shutdown_date = None
        self.shutdown_time = None

        self.load_list = []
        self.backstop_list = []
        self.load_type_list = []
        self.continuity_file_name = cont_file_name
        self.NLET_tracking_file_path = NLET_tracking_file_path

        # List of WSPOW commands
        self.power_command_list = [ 'WSPOW00000', 'WSPOW0002A', 'WSVIDALLDN']


        # Full path to RTS files
        self.RTS = LTCTI_RTS.LTCTI_RTS(os.path.dirname(__file__))

        # List of any events that occurred between a shutdown and
        # The time of first command of the Return to Science load.
        #   - Or really between any 2 times but the above is how it's used.
        self.shutdown_event_list = []

        # Create a Dtype for the Continuity Info array
        self.cont_dtype = [('base_load', '|U20'),
                           ('cont_file', '|U80'),
                           ('load_type', '|U10'),
                           ('load_tofc', '|U25')]

        # Create the SCS-107 command set that's relevant for ACIS
        # Now the SCS-107 on board is an RTS so the times and
        # VCDU's in this data structure are bogus.  Meaningful times
        # will be entered in a copy of this data structure when the commands
        # are integrated into a backstop history
        # NOTE: Never write into this attribute. Always make a copy
        self.scs107_bs_cmds = [
                              # --------------------- SIMTRANS -----------------------------
                               {'cmd': 'SIMTRANS',
                                'date': '1900:001',
                                'msid': None,
                                'params': {'POS': -99616, 'SCS': 108, 'STEP': 1},
                                'paramstr': 'POS= -99616, SCS= 108, STEP= 1',
                                'scs': 108,
                                'step': 1,
                                'time': -1.0,
                                'tlmsid': None,
                                'vcdu': 0000000},

                                # --------------------- AA00000000 -----------------------------
                                {'cmd': 'ACISPKT',               # Stop Science
                                'date': '1900:001',
                                'msid': None,
                                'params': {'CMDS': 3,
                                           'PACKET(40)': 'D80000300030603001300',
                                           'SCS': 107,
                                           'STEP': 1,
                                           'TLMSID': 'AA00000000',
                                           'WORDS': 3},
                                'paramstr': 'TLMSID= AA00000000, CMDS= 3, WORDS= 3, PACKET(40)= D80000300030603001300                   , SCS= 107, STEP= 1',
                                'scs': 107,
                                'step': 1,
                                'time': -1.0,
                                'tlmsid': 'AA00000000',
                                'vcdu': 0000000},

                               # --------------------- AA00000000 -----------------------------
                                {'cmd': 'ACISPKT',               # Stop Science
                                'date': '1900:001',
                                'msid': None,
                                'params': {'CMDS': 3,
                                           'PACKET(40)': 'D80000300030603001300',
                                           'SCS': 107,
                                           'STEP': 2,
                                           'TLMSID': 'AA00000000',
                                           'WORDS': 3},
                                'paramstr': 'TLMSID= AA00000000, CMDS= 3, WORDS= 3, PACKET(40)= D80000300030603001300                   , SCS= 107, STEP= 2',
                                'scs': 107,
                                'step': 2,
                                'time': -1.0,
                                'tlmsid': 'AA00000000',
                                'vcdu': 0000000},
                               # --------------------- WSPOW00000 -----------------------------
                               {'cmd': 'ACISPKT',
                                'date': '1900:001',
                                'msid': None,
                                'params': {'CMDS': 5,
                                           'PACKET(40)': 'D8000070007030500200000000000010000',
                                           'SCS': 107,
                                           'STEP': 3,
                                           'TLMSID': 'WSPOW00000',
                                           'WORDS': 7},
                                'paramstr': 'TLMSID= WSPOW00000, CMDS= 5, WORDS= 7, PACKET(40)= D8000070007030500200000000000010000     , SCS= 107, STEP= 3',
                                'scs': 107,
                                'step': 3,
                                'time': -1.0,
                                'tlmsid': 'WSPOW00000',
                                'vcdu': 0000000} ]

#                                'vcdu': 0000000},


        # --------------------------------------------------------------
        #
        #  Individual Power Commands - Commands issued by themselves
        #
        # --------------------------------------------------------------
        # --------------------- WSPOW00000 -----------------------------
        self.WSPOW00000_bs_cmd = {'cmd': 'ACISPKT',
                                'date': '1900:001',
                                'msid': None,
                                'params': {'CMDS': 5,
                                           'PACKET(40)': 'D8000070007030500200000000000010000',
                                           'SCS': 135,
                                           'STEP': 1,
                                           'TLMSID': 'WSPOW00000',
                                           'WORDS': 7},
                                'paramstr': 'TLMSID= WSPOW00000, CMDS= 5, WORDS= 7, PACKET(40)= D8000070007030500200000000000010000     , SCS= 135, STEP= 1',
                                'scs': 135,
                                'step': 1,
                                'time': -1.0,
                                'tlmsid': 'WSPOW00000',
                                'vcdu': 0000000}

        # --------------------- WSPOW0002A -----------------------------
        self.WSPOW0002A_bs_cmd = {'cmd': 'ACISPKT',
                                'date': '1900:001',
                                'msid': None,
                                'params': {'CMDS': 5,
                                           'PACKET(40)': 'D80000700073E800020000000000001002A',
                                           'SCS': 135,
                                           'STEP': 1,
                                           'TLMSID': 'WSPOW0002A',
                                           'WORDS': 7},
                                'paramstr': 'TLMSID= WSPOW0002A, CMDS= 5, WORDS= 7, PACKET(40)= D80000700073E800020000000000001002A     , SCS= 135, STEP= 1',
                                'scs': 135,
                                'step': 1,
                                'time': -1.0,
                                'tlmsid': 'WSPOW0002A',
                                'vcdu': 0000000}

        # --------------------- WSVIDALLDN -----------------------------
        self.WSVIDALLDN_bs_cmd = {'cmd': 'ACISPKT',
                                'date': '1900:001',
                                'msid': None,
                                'params': {'CMDS': 4,
                                           'PACKET(40)': 'D800005000506050020000000000',
                                           'SCS': 135,
                                           'STEP': 1,
                                           'TLMSID': 'WSVIDALLDN',
                                           'WORDS': 5},
                                'paramstr': 'TLMSID= WSVIDALLDN, CMDS= 4, WORDS= 5, PACKET(40)= D800005000506050020000000000     , SCS= 135, STEP= 1',
                                'scs': 135,
                                'step': 1,
                                'time': -1.0,
                                'tlmsid': 'WSVIDALLDN',
                                'vcdu': 0000000}

        # Create the MP_TARGQUAT command that specifies a maneuver's final
        # position.
        # VCDU's, SCS, and STEP values in this data structure are bogus.  Meaningful times
        # will be entered in a copy of this data structure when the command
        # is integrated into a backstop history
        # NOTE: Never write into this attribute. Always make a copy
                              # --------------------- MP_TARGQUAT ------------------------
        self.MAN_bs_cmds =  {'cmd': 'MP_TARGQUAT',
                                'date': '1900:001',
                                'msid': None,
                                'params': {'CMDS': 8,
                                           'Q1': 1000.0,
                                           'Q2': 2000.0,
                                           'Q3': 3000.0,
                                           'Q4': 4000.0,
                                           'SCS': 135,
                                           'STEP': 1,
                                           'TLMSID': 'AOUPTARQ'},
                                'paramstr': 'POS= -99616, SCS= 108, STEP= 1',
                                'scs': 135,
                                'step': 1,
                                'time': -1.0,
                                'tlmsid': 'AOUPTARQ',
                                'vcdu': 0000000}

        # Create the AOMANUVR command that initiates a maneuver
        # VCDU's, SCS, and STEP values in this data structure are bogus.  Meaningful times
        # will be entered in a copy of this data structure when the command
        # is integrated into a backstop history
        # NOTE: Never write into this attribute. Always make a copy
                              # --------------------- AOMANUVR ------------------------
        self.AOMANUVR_bs_cmd =  {'cmd': 'COMMAND_SW',
                                'date': '1900:001',
                                'msid': 'AOMANUVR',
                                'params': {'HEX': 8034101,
                                           'MSID': 'AOMANUVR',
                                           'SCS': 135,
                                           'STEP': 2,
                                           'TLMSID': 'AOMANUVR'},
                                'paramstr': 'TLMSID= AOMANUVR, HEX= 8034101, MSID= AOMANUVR, SCS= 135, STEP= 1',
                                'scs': 135,
                                'step': 2,
                                'time': -1.0,
                                'tlmsid': 'AOMANUVR',
                                'vcdu': 0000000}

        # Create the AONSMSAF command that initiates a maneuver
        # VCDU's, SCS, and STEP values in this data structure are bogus.  Meaningful times
        # will be entered in a copy of this data structure when the command
        # is integrated into a backstop history
        # NOTE: Never write into this attribute. Always make a copy
                              # --------------------- AOMANUVR ------------------------
        self.AONSMSAF_bs_cmd =  {'cmd': 'COMMAND_SW',
                                'date': '1900:001',
                                'msid': 'AONSMSAF',
                                'params': {'HEX': 9999999,
                                           'MSID': 'AONSMSAF',
                                           'SCS': 135,
                                           'STEP': 2,
                                           'TLMSID': 'AONSMSAF'},
                                'paramstr': 'TLMSID= AONSMSAF, HEX= 9999999, MSID= AONSMSAF, SCS= 135, STEP= 1',
                                'scs': 135,
                                'step': 2,
                                'time': -1.0,
                                'tlmsid': 'AONSMSAF',
                                'vcdu': 0000000}


#-------------------------------------------------------------------------------
#
# method set_backstop_lists - Given the name of a weekly load (e.g. MAR2717A)
#                             and the name of the backstop file for that load,
#                             insert the load name into the beginning of the class
#                             attribute list:  self.load_list and the backstop file
#                             name at the beginning of the self.backstop_list
#
#-------------------------------------------------------------------------------
    def set_backstop_lists(self, load_week, backstop_name, load_type):
        """
        Given the name of a weekly load (e.g. MAR2717A)
        and the name of the backstop file for that load,
        insert the load name into the beginning of the class
        attribute list:  self.load_list and the backstop file
        name at the beginning of the self.backstop_list

        The idea here is to maintain a history of the files and directories
        you used to create the set of backstop commands.

        After being used once, the lists should be cleared by calling
        self.clear_backstop_lists if you are running two or more histories
        in one program
        """
        if load_week is not None:
            self.load_list.insert(0, load_week)

        if backstop_name is not None:
            self.backstop_list.insert(0, backstop_name)

        if load_type is not None:
            self.load_type_list.insert(0, load_type)

#-------------------------------------------------------------------------------
#
# method clear_backstop_lists - Clear out the load and backstop file history lists
#
#-------------------------------------------------------------------------------
    def clear_backstop_lists(self):
        """
        Clearing out  self.load_list and self.backstop_list
        """
        del self.load_list[:]
        del self.backstop_list[:]
        del self.load_type_list[:]

#-------------------------------------------------------------------------------
#
# method print_backstop_lists - Print the load week, load type and backstop file
#                               name  history lists
#
#-------------------------------------------------------------------------------
    def print_backstop_lists(self):
        """
        Print out self.load_list and self.backstop_list
        """
        print(self.load_list)
        print(self.backstop_list)
        print(self.load_type_list)

#-------------------------------------------------------------------------------
#
#  get_backstop_continuity_path
#
#-------------------------------------------------------------------------------
    def get_continuity_file_info(self, oflsdir):
        """ Given an ofls directory, open the Continuity text file
            within the OFLS directory; read the continuity file path,
            from the first line, the type of load the *REVIEW* load is (not
            the continuity load), on the second line.  If the review load
            is not the normal weekly load (e.g. scs-107 or TOO) grab the time
            of interrupt also from the second line.


               - The load in oflsdir was built using another load as Continuity
                 assuming that the build was either a normal weekly load or a TOO
                 load. You wouldn't be using this function if this is a return to
                 science load after a shutdown of some sort.

             Tack the continuity file name (e.g. ACIS-Continuity.txt) to the path. Open up the
             Continuity text file and read the continuity load path.

             Return the  path to that continuity file ofls directory.

             INPUTS: The OFLS "review" directory

            OUTPUTS: The full path to the continuity ofls backstop file
                     The type of load the REVIEW load is.
                     The time of interrupt if the REVIEW load is an interrupt
                     load.

        """
        # `oflsdir` is like <root>/2018/MAY2118/ofls
        oflsdir = Path(oflsdir)

        # Require that oflsdir starts with /data/acis unless the environment
        # variable ALLOW_NONSTANDARD_OFLSDIR is set.
        allow_nonstandard_oflsdir = 'ALLOW_NONSTANDARD_OFLSDIR' in os.environ
        if (not allow_nonstandard_oflsdir
                and oflsdir.parts[:3] != ('/', 'data', 'acis')):
            raise ValueError('--backstop_file must start with /data/acis. To remove this '
                             'restriction set the environment variable '
                             'ALLOW_NONSTANDARD_OFLSDIR to any value.')

        oflsdir_root = oflsdir.parents[2]  # gives <root>

        # Does a Continuity file exist for the input path
        ofls_cont_fn = oflsdir / self.continuity_file_name

        if ofls_cont_fn.exists():

            # Open the Continuity text file in the ofls directory and read the name of the
            # continuity load date (e.g. FEB2017).  then close the file
            ofls_cont_file = open(ofls_cont_fn, 'r')

            # Read the first line...the path to the continuity load. The
            # continuity path in the file is hardwired to a /data/acis path,
            # independent of user-specified `oflsdir` (e.g.
            # /data/acis/LoadReviews/2018/MAY2118/ofls), so if a non-standard
            # OFLS dir path is allowed then fix that by putting the last 3 parts
            # (2018/MAY2118/ofls) onto the oflsdir root.
            pth = Path(ofls_cont_file.readline().strip())
            if allow_nonstandard_oflsdir:
                continuity_load_path = str(Path(oflsdir_root, *pth.parts[-3:]))
            else:
                continuity_load_path = str(pth)

            # Read the entire second line - load type and possibly the interrupt time
            type_line = ofls_cont_file.readline()
            # Split the line
            split_type_line = type_line.split()

            # The review load type is always there, and always the first item on the line
            #  so capture it
            review_load_type = split_type_line[0]

            # If the review load type is not "Normal", grab the interrupt time
            # or set the interrupt time to "None"
            if review_load_type.upper() != "NORMAL":
                interrupt_time = split_type_line[1]
            else:
                interrupt_time = None

            # Done with the file...close it
            ofls_cont_file.close()

            # Return the Continuity load path to the caller.
            return continuity_load_path, review_load_type, interrupt_time
        else:
            return None, None, None

#-------------------------------------------------------------------------------
#
#   get_bs_cmds - Get the backstop commands that live in the OFLS directories.
#                 These always start with the characters "CR"
#
#-------------------------------------------------------------------------------
    def get_bs_cmds(self, oflsdir):
        """
        Given the path to an ofls directory, this method will call the "globfile"
        to obtain the name of the backstop file that represents the built load.
        It then calls kadi.commands.get_cmds_from_backstop to obtain the list of commands
        Review and Continuity loads appear in the ....ofls/ subdir and always
        begin with the characters "CR"

        INPUT: oflsdir = Path to the OFLS directory (string)

        OUTPUT   : bs_cmds = A list of the commands within the backstop file
                             in ofls directory that represents the  built load.
                                -  list of dictionary items

                   bs_name - Name of the Backstop file containing the commands
                              - e.g. CR249_2305.backstop
        """
        backstop_file_path = globfile(os.path.join(oflsdir, 'CR*.backstop'))
#        self.logger.info("GET_BS_CMDS - Using backstop file %s" % backstop_file_path)

        # Extract the name of the backstop file from the path
        bs_name = os.path.split(backstop_file_path)[-1]

        # Read the commands located in that backstop file
        bs_cmds = kadi.commands.get_cmds_from_backstop(backstop_file_path)
        bs_cmds = bs_cmds.as_list_of_dict()
#        self.logger.info("GET_BS_CMDS - Found %d backstop commands between %s and %s" % (len(bs_cmds),
#                                                                                         bs_cmds[0]['date'],
#                                                                                         bs_cmds[-1]['date']))

        return bs_cmds, bs_name

#-------------------------------------------------------------------------------
#
#    get_vehicle_only_bs_cmds - Get the backstop commands that live in the
#                               OFLS directories. These always start with the
#                               characters "VR"
#
#-------------------------------------------------------------------------------
    def get_vehicle_only_bs_cmds(self, oflsdir):
        """
        Given the path to an ofls directory, this method will call the "globfile"
        to obtain the name of the backstop file that represents the built load.
        It then calls kadi.commands.get_cmds_from_backstop to obtain the list of commands
        Vehicle_only loads appear in the ....ofls/vehicle/ subdir and always
        begin with the characters "VR"

        INPUT: oflsdir = Path to the OFLS directory (string)

        OUTPUT   : bs_cmds = A list of the ommands within the backstop file
                             in ofls directory that represents the  built load.
                                -  list of dictionary items

        """
        backstop_file_path = globfile(os.path.join(oflsdir, "vehicle", 'VR*.backstop'))
        self.logger.info('GET_BS_CMDS - Using backstop file %s' % backstop_file_path)

        # Extract the name of the backstop file from the path
        bs_name = os.path.split(backstop_file_path)[-1]

        # Read the commands located in that backstop file and convert to
        # list of dict.
        bs_cmds = kadi.commands.get_cmds_from_backstop(backstop_file_path)
        bs_cmds = bs_cmds.as_list_of_dict()

        self.logger.info('GET_VEHICLE_ONLY_CMDS - Found %d backstop commands between %s and %s'
                         % (len(bs_cmds), bs_cmds[0]['date'], bs_cmds[-1]['date']))

        return bs_cmds, bs_name


#-------------------------------------------------------------------------------
#
# CombineNormal - Combine the Continuity backstop commands with the review load
#                 backstop commands, taking any overlap into account.
#
#                 Call this when you are reviewing a normal week-to-week load
#
#-------------------------------------------------------------------------------
    def CombineNormal(self, cont_bs_cmds, rev_bs_cmds):
        """
         Combine the Continuity backstop commands with the review load
         backstop commands, taking any overlap into account.

         NOTE: By "review" load we mean the load being reviewed this week OR
               a combination of one or more continuity loads with the load
               being reviewed this week.  This routine can be called multiple
               times - tacking continuity loads to the start of the "master list"

               NOTE: There can be one or more commands at the start of the
                     Review load whose time stamp is prior to the time stamp
                     of the last command in the continuity load. All commands
                     are executed because review commands are put is a different
                     SCS slot than the continuity commands.  therefore all commands
                     from both backstop files must be preserved and interleaved
                     correctly

                 Call this when you are reviewing a normal week-to-week load

                 INPUTS: Continuity load backstop file commands
                         Review Load Backstop file

                OUTPUTS: Date-sorted Backstop commands of the combined Continuity
                         and Review loads.
        """

        # Combine the continuity command list with the review command list
        newlist = cont_bs_cmds+rev_bs_cmds

        # This is a list of dicts. Sort the list based upon the Chandra Date
        # string located in the key: "date". This will interleave all the
        # commands correctly.
        self.master_list = sorted(newlist, key=lambda k: k['time'])

        # Return the sorted command list to the caller
        return self.master_list



#-------------------------------------------------------------------------------
#
# CombineTOO - Combine the Continuity backstop commands with the review load
#              backstop commands, without overlap.
#
#                 Call this when you are reviewing a TOO load
#
#-------------------------------------------------------------------------------
    def CombineTOO(self, cont_bs_cmds, rev_bs_cmds):
        """
         Combine the Continuity backstop commands with the review load
         backstop commands, without overlap.  The Continuity load is
         cut at the time of the first command in the TOO load. If there
         are commands whose execution time is precisely the cut time, they
         will be kept and executed.  The possible AOACRSTD pre-TOFC command
         in the TOO load is placed appropriately.

         NOTE: By "review" load we mean the load being reviewed this week OR
               a combination of one or more continuity loads with the load
               being reviewed this week.  This routine can be called multiple
               times - tacking continuity loads to the start of the "master list"

                 Call this when you are reviewing a TOO load or if a
                 continuity load was broken by a TOO load.

                 INPUTS: Continuity load backstop file commands
                         Review Load Backstop file

                OUTPUTS: Backstop commands of the combined Continuity and Review
                         loads.
        """
        # Get all the Continuity commands up to and including the time of first command of
        # the Review load
        # NOTE: This will automatically take care of the fact that one or more of the first commands
        # in the new load will come before the end of commands in the continuity load.
        #
        self.master_list = [cmd for cmd in cont_bs_cmds if cmd['time'] <= rev_bs_cmds[0]['time']]

        # Now concatenate the review load taking all the commands
        # to the master list
        newlist = self.master_list + rev_bs_cmds

        self.master_list = sorted(newlist, key=lambda k: k['time'])

        return self.master_list



#-------------------------------------------------------------------------------
#
# Collect_Shutodwn_Events - Find all the events that occur during a shutdown
#                           which take place between the shutdown time and
#                           Time of First Command of the Resumption of Science
#                           load
#
#-------------------------------------------------------------------------------
    def Collect_Shutdown_Events(self, shutdown_date, resumption_of_science_date):
        """
        At present there are three kinds of events that can occur after a Shutdown
        which have thermal consequences for ACIS:

            1) Maneuver to a differfent pitch (OCC commanded)
            2) LTCTI measurement (Long Term CTI)
            3) FEPs On/ FEPs Off (altering the number of FEPS on for thermal reasons

        This method looks into the NLET tracking file and finds all those events which
        occurred between the Shutdown time and the Time of First Command of the Return to
        Science load.

        It will include the actual STOP or SCS-107 entry in the NLET tracking file because
        the start time is the shutdown date.

        Sometimes the time of shutdown and the time of the next event may be the same time
        The shutdown time comes from the OCC.

        Recall that the time of subsequent events are entered by ACIS Ops personnel.

        Example: NSM at 2019:248:16:51:18.00

                 The NSM results in a maneuver to 90 degrees. ACIS Ops enters the
                 maneuver using the NLET GUI. It's likely that the ACIS Ops person
                 will enter the same time of the shutdown forthe time of the maneuver.

                 This method will return both of those events plus the LTCTI and MANEUVER
                 that followed:

                 ['2019:248:16:51:18.000', 'STOP', 'HRC-S,HETG-OUT,LETG-OUT,47912,OORMPDS,CSELFMT2,DISA']
                 ['2019:248:16:51:18.000', 'MAN', '90.53', '263.10', '-0.3509236000', '0.6512416600', '-0.4674259400', '0.4839935900']

                 ['2019:249:00:22:52', 'LTCTI', '1494', '1_4_CTI', '001:00:00:00']

                 ['2019:249:02:02:00', 'MAN', '158.69', '287.30', '-0.5457521700', '0.2760205900', '-0.1740841000', '0.7717913400']

            So when subsequent methods use this collection one has to remember that the first
            entry is the shutdown event itself.

        """
        # Convert the start and stop dates to seconds
        tstart = DateTime(shutdown_date).secs
        tstop = DateTime(resumption_of_science_date).secs

        # Initialize the event list to empty
        event_list = []

        # The Non Load Event Tracking file is an input so that different
        # users of this module can have different NLET files.
        nletfile = open(self.NLET_tracking_file_path, 'r')

        # Get the first line
        nletline = nletfile.readline()

        # Process each line. If it starts with a # then ignore it - it's a
        # comment
        #
        # for as long as you have input lines......
        while nletline:
            # Check to see if it's a comment line
            if nletline[0] != '#':
                # Not a comment; check to see if it's a MAN line
                # and if so, if the time stamp on the line is between
                # tstart and tstop
                #
                # Split the line into tokens
                splitline = nletline.split()

                # If it's not a "GO line; and
                #    the time stamp is between tstart and tstop
                # then capture the information
                #
                # NOTE: The reason we have to check to see if it's
                # a GO line is that there is only one token in that line
                if (splitline[0] != 'GO') and \
                   (DateTime(splitline[0]).secs >= tstart) and \
                   (DateTime(splitline[0]).secs <= tstop):

                    # We have found an event that affects our
                    # review load. Capture the values
                   event_list.append(splitline)


            # Read the next line - or try to
            nletline = nletfile.readline()

        # You've read all the lines. Close the file.
        nletfile.close()

        # Return the list of any events that were found
        return event_list


#-------------------------------------------------------------------------------
#
# Add_an_SCS107 - Add the commands for an SCS-107 to the Master list
#
#-------------------------------------------------------------------------------
    def Add_an_SCS107(self, ):
        """
        This method adds, to the Master List, ONLY those commands directly
        pertinent to an SCS-107.

        For timing purposes the time of the first scs-107 command is the
        safemode shutdown time  plus 1 second. As this is the SCS-107 that
        occurs due to a safing action, no sorting of the master list is
        required.

         Inputs: shutdown_time (float)
        Outputs: None

        The expected value of the shutdown time should be in Chandra seconds. But just in case
        someone sends in a DOY string I'll convert it to seconds no matter what.
        """
        # Now make a copy of the SCS-107 commands and populate the times. Then
        # concatenate the 107 commands to the master list
        scs107_bs_cmds = copy.deepcopy(self.scs107_bs_cmds)

        # The starting time for the first scs107 command will be 1 second after the input
        # shutdown time
        base_time = DateTime(self.shutdown_time).secs + 1

        # populate the date and time slots of each command incrementing the
        # times by 4 seconds
        for eachcmd in scs107_bs_cmds:
            eachcmd['time'] = base_time
            eachcmd['date'] = DateTime(base_time).date

            # Increment the base time by one second
            base_time += 4.0

        # Now concatenate scs107 commands specific to this event to the master list
        # which has been trimmed to include only those commands that executed prior to
        # the shutdown. No need to sort these at this point
        self.master_list += scs107_bs_cmds



#-------------------------------------------------------------------------------
#
# Process_Power_Command - Process any NLET entry where the value in the "event"
#                         column matches items in self.power_command_list
#
#-------------------------------------------------------------------------------
    def Process_Power_Command(self, p_c_event_list):
        """
        This method takes as input a list of strings which fully describe a
        Power Command that took place.  It processes that power command and adds
        the appropriate command to the MASTER LIST.  Then it sorts the master list.

        A sample of a power command entry in the NLET file:
        #-------------------------------------------------------------------------------
        #       Time        Event         CAP num
        #-------------------------------------------------------------------------------
        2019:260:02:20:00    WSPOW0002A   4567

        ......translates into this list:

        [ '2019:260:02:20:00',    'WSPOW0002A',   '4567']

        ....and that's what the input argument looks like.

        In order to get to this method, the value in the event column needs to be in
        self.power_command_list.

        For WSPOW commands, only one command will be inseerted into the Master list.

        The CAP number is ignored by this method - it's there for tracking purposes.
        """
        # Capture the important values in the power command event
        p_c_start_date = p_c_event_list[0]
        # Convert the power command start date to seconds
        p_c_start_time = DateTime(p_c_start_date).secs
        # Now grab the actual event
        p_c_event = p_c_event_list[1]

        # Now insert the appropriate command
        # Grab the appropriate power command
        if p_c_event == 'WSPOW00000':
            new_power_command = copy.deepcopy(self.WSPOW00000_bs_cmd)
        elif p_c_event == 'WSPOW0002A':
            new_power_command = copy.deepcopy(self.WSPOW0002A_bs_cmd)
        elif p_c_event == 'WSVIDALLDN':
            new_power_command = copy.deepcopy(self.WSVIDALLDN_bs_cmd)

        # Set the date
        new_power_command['date'] =  p_c_start_date

        # Append the Power Command to the Master List
        self.master_list.append(new_power_command)

        # Sort the Master list
        self.master_list = sorted(self.master_list, key=lambda k: k['time'])

#-------------------------------------------------------------------------------
#
# Process_MANEUVER - Process the maneuver described by man_event and add the
#                    appropriate commands to the MASTER LIST
#
#-------------------------------------------------------------------------------
    def Process_MANEUVER(self, man_event):
        """
        This method takes as input a list of strings which fully describe the
        maneuver that took place.  It processes that maneuver and adds
        the appropriate commands to the MASTER LIST.  Then it sorts the master list.

        Inputs: man_list - list extracted from the NLET tracking file which
                           provides all the details necessary to add the
                           appropriate commands to self.master_list

        Output:  Updated Master List

        NOTE: It is assumed that an SCS-107 set of commands was added to the
              master list prior to the call to this routine. Whenever you have
              a full stop, the code must insert the commands for an SCS-107.
              Since there can be more than one maneuver after a shutdown and before
              the Return To Science (RTS) load, the routine that calls this
              method must attend to that chore.

              We do not want an SCS-107 to be added for ALL maneuver events that
              happen during a shutdown.

        """
        # Capture the important values in the event
        man_start_date = man_event[0]
        # Convert the maneuver start date to seconds
        man_start_time = DateTime(man_start_date).secs
        pitch = man_event[2]
        nom_roll = man_event[3]
        # Capture the Quaternions
        q1 =  man_event[4]
        q2 =  man_event[5]
        q3 =  man_event[6]
        q4 =  man_event[7]

            # If this is a legal maneuver, process it
            if pitch != 0.0:
                # Now form and add the command (in SKA.Parse format - i.e. dict) which
                # specifies the MP_TARGQUAT
                new_maneuver = copy.deepcopy(self.MAN_bs_cmds)

                # Set the dates, times and Q's
            new_maneuver['date'] = man_start_date
                new_maneuver['params']['Q1'] = float(q1)
                new_maneuver['params']['Q2'] = float(q2)
                new_maneuver['params']['Q3'] = float(q3)
                new_maneuver['params']['Q4'] = float(q4)
            new_maneuver['time'] = DateTime(man_start_date).secs
                paramstr = 'TLMSID= AOUPTARQ, CMDS= 8, Q1= %s, Q2= %s, Q3= %s, Q4= %s, SCS= 1, STEP= 1' % (q1, q2, q3, q4)
                new_maneuver['paramstr'] = paramstr

                # Tack the maneuver to the Master List
                self.master_list.append(new_maneuver)

            else: # It's a bogus maneuver entry - the user didn't specify good Q's
                self.logger.warning("Bogus Maneuver Entry! Quaternions badly specified: \n"
                                    "Bad Q's: %g, %g, %g, %g " % (q1, q2, q3, q4) +
                                    "...therefore bogus pitch and roll: %g, %g" % (pitch, roll))

            # Create the AOMANUVR command and add it to the Master List.
            # This command actually kicks the maneuver off.
            # NOTE: This is done whether or not the maneuver data was good.
            #       It allows the loop to continue looking for maneuvers
            aoman = copy.deepcopy(self.AOMANUVR_bs_cmd)
        aoman_time = DateTime(man_start_date).secs + 1
            aoman_date = DateTime(aoman_time).date
            aoman['date'] = aoman_date
            aoman['time'] = aoman_time
            # Tack the maneuver to the Master List
            # Tacking the command to the Master list doesn't really do much if there
            # is no AOUPTARQ command. But it allows you to search for subsequent maneuver
            # commands
            self.master_list.append(aoman)

        # Sort the Master list
        self.master_list = sorted(self.master_list, key=lambda k: k['time'])


#-------------------------------------------------------------------------------
#
# Process_LTCTI - Process the ltcti described by ltcti_list and add the
#                    appropriate commands to the MASTER LIST
#
#-------------------------------------------------------------------------------
    def Process_LTCTI(self, ltcti_event):
        """
        This method takes as input a list of strings which fully describe the
        ltcti that took place.  It processes that ltcti and adds
        the appropriate commands to the MASTER LIST. Then it sorts the master list.

        Inputs: man_list - list extracted from the NLET tracking file which
                           provides all the details necessary to add the
                           appropriate commands to self.master_list

        Output:  Updated Master List

        NOTE: It is assumed that an SCS-107 set of commands was added to the
              master list prior to the call to this routine. Whenever you have
              a full stop, the code must insert the commands for an SCS-107.
              Since there can be more than one ltcti after a shutdown and before
              the Return To Science (RTS) load, the routine that calls this
              method must attend to that chore.

              We do not want an SCS-107 to be added for ALL ltcti events that
              happen during a shutdown.

        """
        # Capture the LTCTI information from the event
        RTS_start_date = ltcti_event[0]
        self.RTS.RTS_name = ltcti_event[3]
        self.RTS.CAP_num = ltcti_event[2]
        self.RTS.NUM_HOURS = ltcti_event[4]

        # If an LTCTI run was found, add it to the master list
        if RTS_start_date is not None:
            # Process the specified RTS file and get a time-stamped numpy array of the data
            cmd_list = self.RTS.processRTS(self.RTS.RTS_name, self.RTS.SCS_NUM, self.RTS.NUM_HOURS, RTS_start_date)

            # Now convert the numpy array into SKA.Parse command format which is a list of dicts
            LTCTI_bs_cmds = self.RTS.convert_ACIS_RTS_to_ska_parse(cmd_list)

            # The LTCTI either ran to completion or was interrupted by the Return to Science
            # load.  Trim any LTCTI CLD commands that occurred ON or AFTER the
            # Return to Science Time of First Command.
            trimmed_LTCTI_bs_cmds = self.Trim_bs_cmds_After_Date(self.Date_OFC, LTCTI_bs_cmds)

            # Concatenate the LTCTI commands to the Master list
            self.master_list += trimmed_LTCTI_bs_cmds

            # Sort the Master List
            self.master_list = sorted(self.master_list, key=lambda k: k['time'])


#-------------------------------------------------------------------------------
#
# CombineSTOP - Combine the Continuity backstop commands with the review load
#               backstop commands, when both the Vehicle and Science loads have
#               been stopped.
#
#-------------------------------------------------------------------------------
    def CombineSTOP(self, cont_bs_cmds, rev_bs_cmds, shutdown_date):
        """
         Combine the INTERRUPTED Continuity backstop commands with the review load
         backstop commands, when both the Vehicle and Science loads have
         been stopped.

         The combination is clean - there will be no interleaved ACA
         commands or any other overlap.

         The Continuity load is cut at the time of shutdown. That value is
         found in the Non-Load Event Tracking file.

         The Time of First Command is obtained by the Review load backstop
         file itself.

         There is usually a gap between those two values.

         The Non-Load Event Tracking File is checked for the existence of
         various events between the shutdown time and the Time of First
         Command of the review load.

         These events include:   Maneuvers (OCC commanded maneuvers)
                                 LTCTI (Long Term CTI measurement)
                                 FEPs on or off

         If an entry exists, the relevant event commands are translated into
         backstop commands which are then translated into the Ska format. Then
         they are added to the stream of commands between the Review and Continuity loads.
         If the Time of First Command of the Review load is before the end of the
         LTCTI, then the LTCTI commands which did not occur are chopped.


         NOTE: By "review" load we mean the load being reviewed this week OR
               a combination of one or more continuity loads with the load
               being reviewed this week.  This routine can be called multiple
               times - tacking continuity loads to the start of the "master list"


                 INPUTS: Continuity load backstop file commands
                         Review Load Backstop commands
                         Time of Shutdown


                OUTPUTS: Backstop commands of the combined Continuity and Review
                         loads.
        """
        # Convert shutdown date to seconds:
        self.shutdown_date = shutdown_date
        self.shutdown_time = DateTime(shutdown_date).secs

        # Capture the Time Of First Command (OFC) from the rev_bs_cmds
        self.Date_OFC = rev_bs_cmds[0]['date']
        self.Time_OFC = rev_bs_cmds[0]['time']

        # Trim the Continuity commands list to include only those
        # commands whose excution time is before the shutdown time
        self.master_list = [cmd for cmd in cont_bs_cmds if (cmd['time'] < self.shutdown_time)   ]

        #
        # IMPORTANT: At this point, self.master_list should consist ONLY of
        #            the TRIMMED continuity load

        # Given that we are here due to a Full Stop, we have to insert a set of
        # SCS-107 commands.  Thermal consequences are that ACIS is shut down
        # and not generating heat.

        # Now we need to process any events that appear in the NonLoadEventTracker.txt
        # file whose times are after the stop time, but before the subsequent
        # start science of the review load. These events include:
        #
        #           NSM - pitch change to 90 degrees: ALL Stop
        #           BSH - stuck at some pitch - ALL Stop
        #           MAN - OCC Pitch Maneuver - Move to a new pitch.
        #
        # Next step is to insert SCS-107 commanding into the master list since
        # whatever caused this FULL STOP caused an SCS-107 run

        self.Add_an_SCS107()

        # Collect any and all events that occurred between the Shutdown time and the
        # Time of First Command of the Review load
        collected_events = self.Collect_Shutdown_Events(self.shutdown_date, rev_bs_cmds[0]['time'])

        # Now we have collected any and all events that occurred after the shutdownincluding the shutdown itself.
        # Process each event, tacking the appropriate command lists onto the Master List
        # Remember that the first event in the list is the STOP shutdown itself
        for man_num, eachevent in enumerate(collected_events[1:]):
            # MANEUVER CHECK - Check to see if the event is a MANERUVER
            if eachevent[1] == 'MAN':
                self.Process_MANEUVER(eachevent)
            elif eachevent[1] == 'LTCTI':
                self.Process_LTCTI(eachevent)
            elif eachevent[1] in self.power_command_list:
                self.Process_Power_Command(eachevent)

        # Finally,  concatenate the REVIEW load taking all the commands
        # to the master list.
        # NOTE: In subsequent calls to this method the TOFC of rev_bs_cmds will be
        #       the start of the *assembled* load history. Only in the case of the first
        #       call is the TOFC of rev_bs_cmds also the TOFC of the actual Review Load.
        newlist =  self.master_list + rev_bs_cmds

        # Final sort
        self.master_list = sorted(newlist, key=lambda k: k['time'])

        return self.master_list


#-------------------------------------------------------------------------------
#
# Combine107 - Combine the Continuity backstop commands with the review load
#              backstop commands, without overlap and without including any
#              command after the specified SCS-107 time
#
#                 Call this when you are reviewing an SCS-107 load
#
#-------------------------------------------------------------------------------
    def Combine107(self, cont_bs_cmds, vo_bs_cmds, rev_bs_cmds, shutdown_date):
        """
         Combine the Continuity backstop commands with the review load
         backstop commands, without overlap IN THE SIENCE LOAD and without including any
         SCIENCE command after the specified SCS-107 time.  The Continuity load is
         cut at the time of shutdown.

         THe SCS 107 commands are inserted after the time of shutdown

         LTCTI is looked for in the NLET file and if it exists, it is inserted
         into the history

         There is a check for the  AOACRSTD command in the TOO load and it's
         placed accordingly via a sort.

         Steps:

             1) Trim all commands from the continuity load on and after the
                SCS-107 interrupt time.

             2) Add in the SCS-107 Shutdown commands (SIMTRANS, AA, AA, WSPOW)

             3) Check for a LTCTI

             4) Trim the beginning of the VO load to eliminate that part of the
                continuity load you are keeping.

             5) Trim the back end of the VR file to the Review load start time

             6)  Concatenate the VO remnant to the assembled continuity.

             7) Concatenate the Review Load

         NOTE: By "review" load we mean the load being reviewed this week OR
               a combination of one or more continuity loads with the load
               being reviewed this week.  This routine can be called multiple
               times - tacking continuity loads to the start of the "master list"

                 Call this when you are reviewing a normal week-to-week load

                 INPUTS: Continuity load backstop file commands
                         Review Load Backstop file
                         Time of Shutdown

                OUTPUTS: Backstop commands of the combined Continuity and Review
                         loads.
        """
        # Convert shutdown date to seconds:
        self.shutdown_date = shutdown_date
        self.shutdown_time = DateTime(shutdown_date).secs

        # Capture the Time of First Command from the rev_bs_cmds
        self.Date_OFC = rev_bs_cmds[0]['date']
        self.Time_OFC = rev_bs_cmds[0]['time']

        # STEP 1

        # Trim the Continuity commands list to include only those commands whose excution
        # time is less than the time of first command of the Review load (or assembled bs list)
        self.master_list = [cmd for cmd in cont_bs_cmds if (cmd['time'] < self.shutdown_time)   ]

        # Capture the end of the master list (trimmed continuity )
        # so that you can use it later for the Verhicle Only Load Cut. You want it to
        # be right after the SCS-107 time. But we are adding in any LTCTI runs first.
        # That would move the end time of master list too far. The same thing would happen
        # to the LTCTI if we swapped them in time.
        vo_cut_time = self.master_list[-1]['time']
        vo_cut_date = self.master_list[-1]['date']

        # STEP 2 SCS-107 SIMTRANS AND STOP
        self.Add_an_SCS107()


        # STEP 3 - Collect Events
        # Collect any and all events that occurred between the Shutdown time and the
        # Time of First Command of the Review load
        collected_events = self.Collect_Shutdown_Events(self.shutdown_date, rev_bs_cmds[0]['time'])

        # Now we have collected any and all events that occurred after the shutdownincluding the shutdown itself.
        # Process each event, tacking the appropriate command lists onto the Master List
        # Remember that the first event in the list is the STOP shutdown itself
        for man_num, eachevent in enumerate(collected_events[1:]):
            # MANEUVER CHECK - Check to see if the event is a MANERUVER
            if eachevent[1] == 'MAN':
                self.Process_MANEUVER(eachevent)
            elif eachevent[1] == 'LTCTI':
                self.Process_LTCTI(eachevent)

        # STEP 4
        vo_bs_cmds_trimmed = self.Trim_bs_cmds_Before_Date(vo_cut_date, vo_bs_cmds)

        # STEP 5

        # Trim the VO list by removing all commands that occur on and after the RTS load
        vo_bs_cmds_trimmed = self.Trim_bs_cmds_After_Date(rev_bs_cmds[0]['time'], vo_bs_cmds_trimmed)

        # STEP 6

        # Concatenate the VO list
        self.master_list += vo_bs_cmds_trimmed

        # STEP 7
        # Finally,  concatenate the review load taking all the commands
        # to the master list
        newlist =  self.master_list + rev_bs_cmds

        # sort them
        self.master_list = sorted(newlist, key=lambda k: k['time'])

        # Return the expanded Master List to the caller
        return self.master_list


    #-------------------------------------------------------------------------------
    #
    # WriteCombinedCommands_Debug - Write the combined commands to a file
    #
    #-------------------------------------------------------------------------------
    """
    This method will write the command list out into a file whose path is specified
    in outfile_path. This is a special version of writing commands which is used for
    debug only.  You can either write a new file out, or append TO an existing file.
    This is controlled by the "new_or_append" input argument.
    This method also writes a dividing line into the file to help you identify chunks
    of commands that the combiner code has inserted. This line is written first prior
    to the associated commands.

    So the dividing line is most useful when appending commands

    Don't forget that the cmd_list does not have to be the entire master list of commands.
    You can call this method after you have assembled a spcific command set - like the set
    of commands from an SCS-107, or a LTCTI RTS - and feed that specific command set to a
    file with the "append" argument.

    THEREFORE the resultant file is NOT USABLE in a thermal model.  It's there only
    for you to debug the combiner code appearin in this class.

        INPUTS: command list
                full output file specification
                new_or_append - "new" or "append"
                event_id - a string to be inserted in the middle of the dividing line
                         - the string can be whatever helps you ID the command chunk.

       OUTPUTS: Nothing returned; file written.
    """
    def WriteCombinedCommands_Debug(self, cmd_list, outfile_path, new_or_append = 'new', event_id = 'None'):
        # Write a new file or append to an existing file
        if new_or_append == 'new':
            combofile = open(outfile_path, 'w')
        else:
            combofile = open(outfile_path, 'a')

        # Write out the dividing line
        combofile.write('------------------------------ '+ event_id + ' -----------------------------------\n')

        if event_id != 'MANEUVER':
            # Now put out the set of commands supplied by the caller
            for eachcmd in cmd_list:
                if eachcmd['cmd'] != 'GET_PITCH':
                    cmd_line = eachcmd['date'] + " | %s | %s | %s\n" % (str(eachcmd['vcdu']).zfill(7),
                                                                        eachcmd['cmd'],
                                                                        eachcmd['paramstr'])
                    combofile.write(cmd_line)
        else: # Else it's a single dict entry
            cmd_line = cmd_list['date'] + " | %s | %s | %s\n" % (str(cmd_list['vcdu']).zfill(7),
                                                                        cmd_list['cmd'],
                                                                        cmd_list['paramstr'])
            combofile.write(cmd_line)

        # Done with the file....close it
        combofile.close()

    #-------------------------------------------------------------------------------
    #
    # WriteCombinedCommands - Write the combined commands to a file
    #
    #-------------------------------------------------------------------------------
    """
    This method will write the command list out into a file whose path is specified
    in outfile_path.  Whether or not this is an original command list or a combined
    one it immaterial.

        INPUTS: command list
                full output file specification

       OUTPUTS: Nothing returned; file written.
    """
    def WriteCombinedCommands(self, cmd_list, outfile_path):
        combofile = open(outfile_path, "w")
        for eachcmd in cmd_list:
            if eachcmd['cmd'] != 'GET_PITCH':
                cmd_line = eachcmd['date'] + " | %s | %s | %s\n" % (str(eachcmd['vcdu']).zfill(7),
                                                                    eachcmd['cmd'],
                                                                    eachcmd['paramstr'])
                combofile.write(cmd_line)

        combofile.close()


    #-------------------------------------------------------------------------------
    #
    # Trim_bs_commands_After_Date - Given a list of backstop commands, remove any command
    #                               prior to the specified time
    #
    #                 INPUTS: Chandra date or time (either is ok)
    #                         List of backstop commands.
    #
    #                OUTPUTS: Trimmed list of backstop commands
    #-------------------------------------------------------------------------------
    def Trim_bs_cmds_After_Date(self, trim_date, cmd_list):
        """
        Given a list of backstop commands, remove any command ON or *AFTER* the specified time

        This can be confusing - you want to throw away all those commands that occur
        ON or AFTER the specified time. You want to keep everything that occurs BEFORE the
        cutoff time.

        INPUTS: Chandra date or time (either is ok)

                List of backstop commands.

        OUTPUTS: Trimmed list of backstop commands

        NOTE: John Scott assures me that they try to schedule the first
              command of the new load to be either between loads or in between two
              commands in the old load.  Also, they will not lengthen an observation
              by having the stop science come much later than it did in the broken load.
              They tried that once and there were myriad complaints.
        """
        # Convert the date into seconds. If it was already seconds this will not hurt.
        trim_time = DateTime(trim_date).secs

        # Now trim the list by including only those backstop commands that occur after the time.
        trimmed_list = [cmd for cmd in cmd_list if cmd['time'] < trim_time]

        # Return the trimmed list
        return trimmed_list



    #-------------------------------------------------------------------------------
    #
    # Trim_bs_cmds_Before_Date - Given a list of backstop commands, remove any command
    #                            prior to the specified time
    #
    #                 INPUTS: Chandra date or time (either is ok)
    #                         List of backstop commands.
    #
    #                OUTPUTS: Trimmed list of backstop commands
    #-------------------------------------------------------------------------------
    def Trim_bs_cmds_Before_Date(self, trim_date, cmd_list):
        """
        Given a list of backstop commands, remove any command *PRIOR* to the specified time

        This can be confusing - you want to throw away all those commands that occur BEFORE
        the specified time. You want to keep everything that occurs ON or AFTER the
        cutoff time.

        INPUTS: Chandra date or time (either is ok)
                List of backstop commands.

        OUTPUTS: Trimmed list of backstop commands

        NOTE: John Scott assures me that you want to do a "less than" comparison
              on times and not less than or equal. Also they try to schedule the first
              command of the new load to be either between loads or in between two
              commands in the old load.  Also, they will not lengthen an observation
              by having the stop science come much later than it did in the broken load.
              They tried that once and there were myriad complaints.
        """
        # Convert the date into seconds. If it was already seconds this will not hurt.
        trim_time = DateTime(trim_date).secs

        # Now trim the list by KEEPING only those backstop commands that occur AFTER the time.
        trimmed_list = [cmd for cmd in cmd_list if cmd['time'] >= trim_time]

        # Return the trimmed list
        return trimmed_list


    #-------------------------------------------------------------------------------
    #
    # Backchain - Given a base directory, a starting load (base_load), and a
    #             chain length, this method will successively backtrack through the
    #             Continuity text files of each load starting with the
    #             base load, and record the continuity information in a numpy array
    #
    #                 INPUTS: base_dir  (e.g. '/data/acis/LoadReviews/2017/')
    #                         base_load (e.g. 'AUG3017')
    #                         chain_length - the number of backtracks you want
    #                         to make.
    #
    #                OUTPUTS: Array of records for each back chain through the
    #                         Continuity files.
    #
    #    VITALLY IMPORTANT!!!!!!! The January 30, 2017 load was the FIRST LOAD
    #                             to have the Continuity text file stored.
    #                             Therefore you cannont back Chain further beyond
    #                             The January 30th load.
    #
    #-------------------------------------------------------------------------------
    def BackChain(self, base_load_dir, chain_length):
        """
        Given a full base load directory, and a
        chain length, this method will successively backtrack through the
        Continuity text files of each load starting with the
        base load, and record the continuity information in a numpy array

        INPUTS: base_dir  (e.g. '/data/acis/LoadReviews/2017/AUG3017/ofls')

                              chain_length - the number of backtracks you want
                             to make.

        OUTPUTS: Array of records for each back chain through the
                 Continuity files.

                 Example, for inputs:
                                     base_dir = '/data/acis/LoadReviews/2017/AUG3017/ofls'
                                     chain_length = 4

                                      The output would be:
     ('AUG3017', '/data/acis/LoadReviews/2017/AUG2817/oflsb', 'TOO', '2017:242:23:35:01.28')
     ('AUG2817', '/data/acis/LoadReviews/2017/AUG2517/oflsc', 'Normal', 'None')
     ('AUG2517', '/data/acis/LoadReviews/2017/AUG1917/oflsa', 'TOO', '2017:237:03:30:01.28')
     ('AUG1917', '/data/acis/LoadReviews/2017/AUG1417/oflsb', 'TOO', '2017:231:15:42:18.91')

        VITALLY IMPORTANT!!!!!!! The January 30, 2017 load was the FIRST LOAD
                                 to have the Continuity text file stored.
                                 Therefore you cannont back Chain further beyond
                                 The January 30th load.

        """
        # Create an empty load chain array
        load_chain = np.array([], dtype=self.cont_dtype)

        # Extract the base load week from the full path
        # This is entered in the first column of the resultant array
        base_load_week = os.path.split(base_load_dir)[-2]
        # Extract the continuity info from the Load week.
        # REMEMBER: this information is with regard to the
        # input load. It's the Continuity file that leads to the
        # Review/input; whether the review/input load is Normal, TOO
        # or SCS-107; and if the latter two - what the Time of
        # First Command is.
        continuity_info = self.get_continuity_file_info(base_load_dir)

        # Count the number of links you have added because the
        # continuity file exists. Use that to know when you've collected enough
        chain_links = 0

        # What we want to do is keep tacking continuity info onto the array
        # until we get the number of entries asked for OR we run into a load
        # week that does not have a continuity file. In the latter case we
        # want to exit gracefully and give the user what we have (if anything).
        #
        # You've done one fetch. It's either loaded with continuity info
        # (if ACIS-Continuity.txt exists in the directory) or 3 None's were
        # returned.
        #
        # If you have not exceeded the requested chain length and
        # there is continuity info, tack it onto the load chain array
        while len(load_chain) < chain_length and continuity_info[0] is not None:

            continuity_load_path = continuity_info[0]

            # Load up an entry in the array for this load
            load_chain = np.r_[load_chain,
                                np.array( [ (base_load_week,
                                          continuity_load_path,
                                          continuity_info[1],
                                          continuity_info[2]) ],
                                dtype = self.cont_dtype)  ]

            # Try to do it again

            # Extract the base load week from the returned full path
            # This is entered in the first column of the resultant array
            base_load_week = continuity_load_path.split('/')[-2]

            # Get the continuity info for this new week
            continuity_info = self.get_continuity_file_info(continuity_load_path)


        # Return the array of back chains
        return load_chain


    #-------------------------------------------------------------------------------
    #
    #
        #
                #
                #
    #-------------------------------------------------------------------------------


    #-------------------------------------------------------------------------------
    #
    # write_back_chain_to_pickle - Given a backchain, write the contents of the
    #                              backchain to a pickle file. If this is ACIS format
    #                              the full path to the OFLS directory is written out.
    #
    #-------------------------------------------------------------------------------
    def write_back_chain_to_pickle(self, format='ACIS',
                                   file_path='/home/gregg/DPAMODEL/dpa_check/Chain.p', chain=None):
        """
        Given a backchain, and a file path, write the contents of the backchain to a pickle file.
        If the "format" argumen is "ACIS' then the full path to the ACIS ofls directory is
        written out. If the format is anything except ACIS, then the load week
        (e.g. JAN0917) is extracted from the second column in the backchain array and only that is
        written out.

        The latter option allows users with their own backstop tarball directory structure to
        utilize the backchain and obtain backstop command files.

        The format of the backchain must be a numpy array whose DTYPE is self.cont_dtype

        Note that the results of Backchain can be an enpty array. If so, the resultant pickle
        file will be empty.
        """
        # If the format is ACIS then write out the array as is - meaning that the
        # 'cont_file' column will specify a full path such as:
        #         '/data/acis/LoadReviews/2017/OCT0917/ofls'
        if format == 'ACIS':
            with open(file_path, 'w') as outfile:
                chain.dump(outfile)
            outfile.close()
        else: # Just write out the continuity week name from the second column (e.g. OCT0917)
            # Create an empty general array of detype self.cont_dtype
            gen_chain = np.array( [], dtype = self.cont_dtype)

            # Run through the chain, modify the second column, and append it to the general
            # array
            for eachload in chain:
                # Extract just the load week from the full continuity path
                continuity_week = eachload[1].split('/')[-2]

                # Append the  entry in the array for this load
                gen_chain = np.r_[gen_chain,
                                  np.array( [ (eachload[0],      # Base load
                                               continuity_week,  # Continuity load name
                                               eachload[2],      # Base load type
                                               eachload[3]) ],   # Interrupt time or None
                                           dtype = self.cont_dtype)  ]
            # Now you've created the generalized load array so dump that
            with open(file_path, 'w') as outfile:
                gen_chain.dump(outfile)
            outfile.close()

    #-------------------------------------------------------------------------------
    #
    # write_back_chain_to_txt - Given a backchain, write the contents of the
    #                           backchain to a text file. If this is ACIS format
    #                           the full path to the OFLS directory is written out.
    #
    #-------------------------------------------------------------------------------
    def write_back_chain_to_txt(self, format='ACIS',
                                file_path='/home/gregg/DPAMODEL/dpa_check/Chain.txt', chain=None):
        """
        Given a backchain, and a file path, write the contents of the backchain to a text file.
        If the "format" argumen is "ACIS' then the full path to the ACIS ofls directory is
        written out. If the format is anything except ACIS, then the load week
        (e.g. JAN0917) is extracted from the second column in the backchain array and only that is
        written out.

        The latter option allows users with their own backstop tarball directory structure to
        utilize the backchain and obtain backstop command files.

        The format of the backchain must be a numpy array whose DTYPE is self.cont_dtype

        Note that the results of Backchain can be an enpty array. If so, the resultant pickle
        file will be empty.
        """
        # If the format is ACIS then write out the array as is - meaning that the
        # 'cont_file' column will specify a full path such as:
        #         '/data/acis/LoadReviews/2017/OCT0917/ofls'
        if format == 'ACIS':
            np.savetxt(file_path, chain, fmt = '%s', newline = '\n')

        else: # Just write out the continuity week name from the second column (e.g. OCT0917)
            # Create an empty general array of detype self.cont_dtype
            gen_chain = np.array( [], dtype = self.cont_dtype)

            # Run through the chain, modify the second column, and append it to the general
            # array
            for eachload in chain:
                # Extract just the load week from the full continuity path
                continuity_week = eachload[1].split('/')[-2]

                # Append the  entry in the array for this load
                gen_chain = np.r_[gen_chain,
                                  np.array( [ (eachload[0],      # Base load
                                               continuity_week,  # Continuity load name
                                               eachload[2],      # Base load type
                                               eachload[3]) ],   # Interrupt time or None
                                           dtype = self.cont_dtype)  ]
            # Now you've created the generalized load array so dump that
            np.savetxt(file_path, gen_chain, fmt = '%s', newline = '\n')

    #-------------------------------------------------------------------------------
    #
    # read_back_chain_from_txt - Given a path to a txt file written by
    #                            write_back_chain_to_txt, read the contents of the
    #                            text file and store it in the specified array
    #                            variable.
    #
    #-------------------------------------------------------------------------------
    def read_back_chain_from_txt(self, file_path='/home/gregg/DPAMODEL/dpa_check/Chain.txt'):
        """
        Given a path to a txt file written by write_back_chain_to_txt,
        read the contents of the text file and store it in an array
        using DTYPE: self.cont_dtype
        """
        # Read the file
        chain = np.loadtxt(file_path, self.cont_dtype)
        # Return the array
        return chain
