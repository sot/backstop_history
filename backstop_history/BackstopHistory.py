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
import copy
import glob
import numpy as np
import os

import LTCTI_RTS

import Ska.ParseCM
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


class BackstopHistory(object):

    def __init__(self, cont_file_name='ACIS-Continuity.txt', 
                 NLET_tracking_file_path='/data/acis/LoadReviews/NonLoadTrackedEvents.txt'):
        self.master_list = []
        self.rev_to_take = []
        self.load_list = []
        self.backstop_list = []
        self.load_type_list = []
        self.continuity_file_name = cont_file_name
        self.NLET_tracking_file_path = NLET_tracking_file_path
        # Full path to RTS files
        self.RTS = LTCTI_RTS.LTCTI_RTS(os.path.dirname(__file__))


        # Create a Dtype for the Continuity Info array
        self.cont_dtype = [('base_load', '|S20'),
                           ('cont_file', '|S80'), 
                           ('load_type', '|S10'), 
                           ('load_tofc', '|S25')]

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
        Print out  self.load_list and self.backstop_list
        """
        print self.load_list[:]
        print self.backstop_list[:]
        print self.load_type_list[:]
    
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
        # Does a Continuity file exist for the input path
        if os.path.isfile(oflsdir+'/'+self.continuity_file_name):
    
            # Open the Continuity text file in the ofls directory and read the name of the 
            # continuity load date (e.g. FEB2017).  then close the file
            ofls_cont_file = open(oflsdir+'/'+self.continuity_file_name, 'r')
            # Read the first line...the pat to the continuity load
            continuity_load_path = ofls_cont_file.readline()[:-1]
    
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
        It then calls Ska.ParseCM.read_backstop to obtain the list of commands
        Review and Continuity loads appear in the ....ofls/ subdir and always
        begin with the characters "CR"

        INPUT: oflsdir = Path to the OFLS directory (string)
    
        OUTPUT   : bs_cmds = A list of the ommands within the backstop file 
                             in ofls directory that represents the  built load.
                                -  list of dictionary items
    
        """
        backstop_file_path = globfile(os.path.join(oflsdir, 'CR*.backstop'))
        print'    GET_BS_CMDS - Using backstop file %s' % backstop_file_path
    
        # Extract the name of the backstop file from the path
        bs_name = backstop_file_path.split('/')[-1]
    
        # Read the commands located in that backstop file
        bs_cmds = Ska.ParseCM.read_backstop(backstop_file_path)
        print '    GET_BS_CMDS - Found %d backstop commands between %s and %s' % (len(bs_cmds), bs_cmds[0]['date'], bs_cmds[
    -1]['date'])
    
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
        It then calls Ska.ParseCM.read_backstop to obtain the list of commands
        Vehicle_only loads appear in the ....ofls/vehicle/ subdir and always
        begin with the characters "VR"

        INPUT: oflsdir = Path to the OFLS directory (string)
    
        OUTPUT   : bs_cmds = A list of the ommands within the backstop file 
                             in ofls directory that represents the  built load.
                                -  list of dictionary items
    
        """
        backstop_file_path = globfile(os.path.join(oflsdir+'/vehicle/', 'VR*.backstop'))
        print'    GET_BS_CMDS - Using backstop file %s' % backstop_file_path
    
        # Extract the name of the backstop file from the path
        bs_name = backstop_file_path.split('/')[-1]
    
        # Read the commands located in that backstop file
        bs_cmds = Ska.ParseCM.read_backstop(backstop_file_path)
        print '    GET_VEHICLE_ONLY_CMDS - Found %d backstop commands between %s and %s' % (len(bs_cmds), bs_cmds[0]['date'], bs_cmds[
    -1]['date'])
    
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
#        self.master_list = sorted (newlist, key=lambda k: k['date'])
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
# CombineSTOP - Combine the Continuity backstop commands with the review load
#               backstop commands, when both the Vehicle and Science loads have
#               been stopped.
#
#-------------------------------------------------------------------------------
    def CombineSTOP(self, cont_bs_cmds, rev_bs_cmds, shutdown_date):
        """
         Combine the Continuity backstop commands with the review load
         backstop commands, when both the Vehicle and Science loads have
         been stopped.

         The combination is clean - there will be no interleaved ACA 
         commands or any other overlap. 

         The Continuity load is cut at the time of shutdown. That value is
         found in the Non-Load Event Tracking file.

         The Time of First Command is obtained by the Review load backstop
         file itself. 

         There is usually a gap between those two values. 

         The Non-Load Event Tracking File is checked for the existence of a LTCTI
         run after the shutdown but before the start of the Review Load. If an
         entry exists, the commands from the relevant CLD file are translated into
         backstop commands which are then translated into the Ska format. Then
         they are added to the stream of commands between the Review and Continuity loads.
         If the Time of First Command of the Review load is before the end of the 
         LTCTI, then the LTCTI commands which did not occur are chopped.


         NOTE: By "review" load we mean the load being reviewed this week OR
               a combination of one or more continuity loads with the load
               being reviewed this week.  This routine can be called multiple
               times - tacking continuity loads to the start of the "master list"

                 
                 INPUTS: Continuity load backstop file commands
                         Review Load Backstop file
                         Review Load Vehicle Only Backstop file
                         Time of Shutdown
                         

                OUTPUTS: Backstop commands of the combined Continuity and Review
                         loads.
        """
        # Convert shutdown date to seconds:
        shutdown_time = DateTime(shutdown_date).secs

        # Capture the Time of First Command from the rev_bs_cmds
        # This is to make the code more self-documenting.
        Date_of_First_Command = rev_bs_cmds[0]['date']
        Time_of_First_Command = rev_bs_cmds[0]['time']

        # Trim the Continuity commands list to include only those 
        # commands whose excution time is before the shutdown time
        self.master_list = [cmd for cmd in cont_bs_cmds if (cmd['time'] < shutdown_time)   ]

        #
        # IMPORTANT: At this point, self.master_list should consist ONLY of
        #            the TRIMMED continuity load

        # Now make a copy of the SCS-107 commands and populate the times. Then
        # concatenate the 107 commands to the master list
        scs107_bs_cmds = copy.deepcopy(self.scs107_bs_cmds)
        
        # The starting time for the first scs107 command will be 1 second after the last
        # command in the TRIMMED master list
        base_time = self.master_list[-1]['time'] + 1

        # populate the date and time slots of each command incrementing the times by one second
        for eachcmd in scs107_bs_cmds:
            eachcmd['time'] = base_time
            eachcmd['date'] = DateTime(base_time).date

            # Increment the base time by one second
            base_time += 1.0
        
        # Now concatenate scs107 commands specific to this event to the master list
        # which has been trimmed to include only those commands that executed prior to 
        # the shutdown. No need to sort these at this point
        self.master_list += scs107_bs_cmds

        # Now we need to process any events that appear in the NonLoadEventTracker.txt
        # file whose times are after the stop time, but before the subsequent
        # start science of the review load. These events include:
        #
        #           NSM - pitch change to 90 degrees: ALL Stop
        #           BSH - stuck at some pitch - ALL Stop
        #     OCC Pitch Maneuver - Move to a new pitch.
        #
        # Next we determine if a MAN (maneuver) entry exists in the NLET file 
        # between the shutdown_date and the time of the first command in "rev_bs_cmds"
        # These could be both NSM and OCC-commanded pitch changes.
        MAN_date, pitch, roll, q1, q2, q3, q4 = self.FindMANs(shutdown_date, rev_bs_cmds[0]['time'])

        while MAN_date is not None:
            # If this is a legal maneuver, process it
            if pitch != 0.0:
                print "\n     MANEUVER FOUND!", MAN_date
                # Now form and add the command (in SKA.Parse format - i.e. dict) which
                # specifies the MP_TARGQUAT 
                new_maneuver = copy.deepcopy(self.MAN_bs_cmds)
    
                # Set the dates, times and Q's
                new_maneuver['date'] = MAN_date
                new_maneuver['params']['Q1'] = float(q1)
                new_maneuver['params']['Q2'] = float(q2)
                new_maneuver['params']['Q3'] = float(q3)
                new_maneuver['params']['Q4'] = float(q4)
                new_maneuver['time'] = DateTime(MAN_date).secs
                paramstr = 'TLMSID= AOUPTARQ, CMDS= 8, Q1= '+str(q1)+', Q2= '+str(q2)+', Q3= '+str(q3)+', Q4= '+str(q4)+'SCS= 1, STEP= 1'
                new_maneuver['paramstr'] = paramstr
    
                # Tack the maneuver to the Master List
                self.master_list += [new_maneuver]
      
            else: # It's a bogus maneuver entry - the user didn't specify good Q's
                print "Bogus Maneuver Entry! Quaternions badly specified: \n"
                print "Bad Q's: ", q1, q2, q3, q4
                print "...therefore bogus pitch and roll: ", pitch, roll

            # Create the AOMANUVR command and add it to the Master List.
            # This command actually kicks the maneuver off.
            # NOTE: This is done whether or not the maneuver data was good.
            #       It allows the loop to continue looking for maneuvers
            aoman = copy.deepcopy(self.AOMANUVR_bs_cmd)
            aoman_time = DateTime(MAN_date).secs + 1
            aoman_date = DateTime(aoman_time).date 
            aoman['date'] = aoman_date
            aoman['time'] = aoman_time
            # Tack the maneuver to the Master List
            # Tacking the command to the Master list doesn't really do much if there
            # is no AOUPTARQ command. But it allows you to search for subsequent maneuver 
            # commands
            self.master_list += [aoman]

            # See if there is another one between the one you found and the beginning of the
            # assembled load
            MAN_date, pitch, roll, q1, q2, q3, q4 = self.FindMANs(aoman['time'], rev_bs_cmds[0]['time'])

        # Next we determine if there is a Long Term CTI run that was done between the
        # last time in the master list and the time of the first command in "rev_bs_cmds"
        RTS_start_date, self.RTS.RTS_name, self.RTS.CAP_num,  self.RTS.NUM_HOURS = self.FindLTCTIrun(shutdown_date , rev_bs_cmds[0]['time'])

        # If an LTCTI run was found, add it to the master list
        if RTS_start_date is not None:
            # Process the specified RTS file and get a time-stamped numpy array of the data
            cmd_list = self.RTS.processRTS(self.RTS.RTS_name, self.RTS.SCS_NUM, self.RTS.NUM_HOURS, RTS_start_date)
    
            # Now convert the numpy array into SKA.Parse command format which is a list of dicts
            LTCTI_bs_cmds = self.RTS.convert_ACIS_RTS_to_ska_parse(cmd_list)
   
            # The LTCTI either ran to completion or was interruopted by the Return to Science
            # load.  Trim any LTCTI CLD commands that occurred ON or AFTER the
            # Return to Science Time of First Command. 
            trimmed_LTCTI_bs_cmds = self.Trim_bs_cmds_After_Date(Date_of_First_Command, LTCTI_bs_cmds)
    
            # Concatenate the LTCTI commands to the Master list
            self.master_list += trimmed_LTCTI_bs_cmds
    
        # Finally,  concatenate the review load taking all the commands
        # to the master list.
        # NOTE: In subsequent calls to this method the TOFC of rev_bs_cmds will be
        #       the start of the *assembled* load history. Only in the case of the first
        #       call is the TOFC of rev_bs_cmds also the TOFC of the actual Review Load.
        newlist =  self.master_list + rev_bs_cmds

        # sort them
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

             7) Concatenate the REview Load

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
        shutdown_time = DateTime(shutdown_date).secs

        # Capture the Time of First Command from the rev_bs_cmds
        # This is to make the code more self-documenting.
        Date_of_First_Command = rev_bs_cmds[0]['date']
        Time_of_First_Command = rev_bs_cmds[0]['time']

        # STEP 1

        # Trim the Continuity commands list to include only those commands whose excution
        # time is less than the time of first command of the Review load (or assembled bs list)
        self.master_list = [cmd for cmd in cont_bs_cmds if (cmd['time'] < shutdown_time)   ]

        # Capture the end of the master list (trimmed continuity )
        # so that you can use it later for the Verhicle Only Load Cut. You want it to
        # be right after the SCS-107 time. But we are adding in any LTCTI runs first.
        # That would move the end time of master list too far. The same thing would happen
        # to the LTCTI if we swapped them in time.
        vo_cut_time = self.master_list[-1]['time']
        vo_cut_date = self.master_list[-1]['date']

        # STEP 2 SCS-107 SIMTRANS AND STOP

        # Now make a copy of the SCS-107 commands and populate the times. These
        # are: SIMTRANS, AA00, AA00, WSPOW
        scs107_bs_cmds = copy.deepcopy(self.scs107_bs_cmds)
        # The starting time for the first scs107 command will be 1 second after the last
        # command in the TRIMMED master list
        base_time = self.master_list[-1]['time'] + 1

        # populate the date and time slots of each command incrementing the times by one second
        for eachcmd in scs107_bs_cmds:
            eachcmd['time'] = base_time
            eachcmd['date'] = DateTime(base_time).date
            # Increment the base time by one second
            base_time += 1.0
        
        # Now concatenate scs107 commands specific to this event to the master list
        # which has been trimmed to include only those commands that executed prior to 
        # the shutdown. No need to sort these at this point, because the base time was
        # the last time of the trimmed continuity load plus 1 second.
        self.master_list += scs107_bs_cmds


        # STEP 3
        print "\n    STEP - 3 Check for a LTCTI between: ",self.master_list[-1]['date'] , rev_bs_cmds[0]['date']
        # Next we determine if there is a Long Term CTI run that was done between the
        # last time in the master list and the time of the first command in "rev_bs_cmds"
        RTS_start_date, self.RTS.RTS_name, self.RTS.CAP_num,  self.RTS.NUM_HOURS = self.FindLTCTIrun(shutdown_date , rev_bs_cmds[0]['time'])

        # If an LTCTI run was found, add it to the master list
        if RTS_start_date is not None:
            # Process the specified RTS file and get a time-stamped numpy array of the data
            cmd_list = self.RTS.processRTS(self.RTS.RTS_name, self.RTS.SCS_NUM, self.RTS.NUM_HOURS, RTS_start_date)
    
            # Now convert the numpy array into SKA.Parse command format which is a list of dicts
            LTCTI_bs_cmds = self.RTS.convert_ACIS_RTS_to_ska_parse(cmd_list)
    
            # The LTCTI either ran to completion or was interruopted by the Return to Science
            # load.  Trim any LTCTI CLD commands that occurred ON or AFTER the
            # Return to Science Time of First Command. 
            trimmed_LTCTI_bs_cmds = self.Trim_bs_cmds_After_Date(Date_of_First_Command, LTCTI_bs_cmds)
    
            # Concatenate the LTCTI run to the master list
            self.master_list += trimmed_LTCTI_bs_cmds
            # Sort the master list
            self.master_list = sorted(self.master_list, key=lambda k: k['time'])
    
        # STEP 4
        vo_bs_cmds_trimmed = self.Trim_bs_cmds_Before_Date(vo_cut_date, vo_bs_cmds)

        # STEP 5

        # Concatenate the VO list
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
    # WriteCombinedCommands - Write the combined commands to a file
    #
    #-------------------------------------------------------------------------------
    """
    This method will write the command list out into a file whose path is specified
    in outfile_path.  Whether or not this is an original command list or a comboned
    one it immaterial.

        INPUTS: command list
                full output file specification

       OUTPUTS: Nothing returned; file written.

    """
    def WriteCombinedCommands(self, cmd_list, outfile_path):
        combofile = open(outfile_path, "w")
        for eachcmd in cmd_list:
            if eachcmd['cmd'] != 'GET_PITCH':
                cmd_line = eachcmd['date']+" | "+str(eachcmd['vcdu']).zfill(7)+" | "+eachcmd['cmd']+" | "+eachcmd['paramstr']+"\n"
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
        load_chain = np.array( [], dtype = self.cont_dtype)

        # Extract the base load week from the full path
        # eliminate a trailing '/' if it's there. 
        # This is entered in the first column of the resultant array
        if base_load_dir[-1] == '/':
            base_load_week = base_load_dir.split('/')[-3]
        else:
            base_load_week = base_load_dir.split('/')[-2]
        
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
    # FindLTCTIrun - Given a path to a Non Load Event Tracking file, a start time
    #                 and a stop time, search the Tracking file for any Long Term 
    #                 CTI run (LTCTI) that occurred between the start and stop times.
    #
    #-------------------------------------------------------------------------------
    def FindLTCTIrun(self, tstart, tstop):
        """
        Given a path to a Non Load Event Tracking file, a start time
        and a stop time, search the Tracking file for any Long Term 
        CTI run (LTCTI) that occurred between the start and stop times.
    
        What you want to use for Start and Stop times are the SCS-107
        times for tstart and the time of first command for the replan load

        The path to the Non Load Event Tracking file (NLET) is a constructor argument
        so that users can have their own version of the file. However the 
        format of the file is fixed and this method expects a certain format.
        """
        # Convert the input tstart and tstop to seconds - this allows the
        # user to input either seconds or DOY format - whichever is more
        # convenient.
        tstart = DateTime(tstart).secs
        tstop = DateTime(tstop).secs
    
        # Initialize the return values to None
        ltcti_date = None
        ltcti_rts_file = None
        ltcti_cap_number = None
        ltcti_duration = None

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
                # Not a comment; check to see if it's a LTCTI line
                # and if so, if the time stamp on the line is between 
                # tstart and tstop
                #
                # Split the line into tokens
                splitline = nletline.split()
                # If it's not a "GO line; is an LTCTI line
                # and the time stamp is between tstart and tstop
                # then capture the information
                #
                # NOTE: The reason we have to check to see if it's
                # a GO line is that ther eis only one token in that line
                if (splitline[0] != 'GO') and \
                   (splitline[1] == 'LTCTI') and \
                   (DateTime(splitline[0]).secs >= tstart) and \
                   (DateTime(splitline[0]).secs <= tstop):
    
                    # We have found a LTCTI event that affects our
                    # review load. Capture the values
                    ltcti_date = splitline[0]
                    ltcti_rts_file = splitline[3]
                    ltcti_cap_number = splitline[2]
                    # NEW
                    ltcti_duration = splitline[4]

            # Read the next line - or try to
            nletline = nletfile.readline()
    
        # You've read all the lines. Close the file.
        nletfile.close()
    
        # Return items from any found netline; or Nones if 
        # no LTCTI line matched the requirements.
        return ltcti_date, ltcti_rts_file, ltcti_cap_number, ltcti_duration




    #-------------------------------------------------------------------------------
    #
    # FindMANs - Given a path to a Non Load Event Tracking file, a start time
    #            and a stop time, search the Tracking file for any maneuver
    #             that occurred between the start and stop times.
    #
    #-------------------------------------------------------------------------------
    def FindMANs(self, dstart, dstop):
        """
        Using the path to a Non Load Event Tracking file, a start time
        and a stop time, search the Tracking file for any maneuver entry
        that occurred between the start and stop times.
    
        Generally what you want to use for Start and Stop times are the SCS-107
        times for tstart and the time of first command for the replan load, or
        the start time for the review plus any accumulated loads or events.

        The path to the Non Load Event Tracking file (NLET) is a constructor argument
        so that users can have their own version of the file. However the 
        format of the file is fixed and this method expects a certain format.
        """
        # Convert the input tstart and tstop to seconds - this allows the
        # user to input either seconds or DOY format - whichever is more
        # convenient.
        tstart = DateTime(dstart).secs
        tstop = DateTime(dstop).secs
    
        # Initialize the return values to None
        MAN_date = None
        MAN_pitch = None
        MAN_roll = None
        MAN_q1 = None
        MAN_q2 = None
        MAN_q3 = None
        MAN_q4 = None
    
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
                # If it's not a "GO line; and is a MAN line
                # and the time stamp is between tstart and tstop
                # then capture the information
                #
                # NOTE: The reason we have to check to see if it's
                # a GO line is that there is only one token in that line
                if (splitline[0] != 'GO') and \
                   (splitline[1] == 'MAN') and \
                   (DateTime(splitline[0]).secs >= tstart) and \
                   (DateTime(splitline[0]).secs <= tstop):
    
                    # We have found a MAN event that affects our
                    # review load. Capture the values
                    MAN_date = splitline[0]
                    MAN_pitch = splitline[2]
                    MAN_roll = splitline[3]
                    MAN_q1 = splitline[4]
                    MAN_q2 = splitline[5]
                    MAN_q3 = splitline[6]
                    MAN_q4 = splitline[7]
     
            # Read the next line - or try to
            nletline = nletfile.readline()
    
        # You've read all the lines. Close the file.
        nletfile.close()
    
        # Return items from any found netline; or Nones if 
        # no LTCTI line matched the requirements.
        return MAN_date, MAN_pitch, MAN_roll, MAN_q1, MAN_q2, MAN_q3, MAN_q4


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
