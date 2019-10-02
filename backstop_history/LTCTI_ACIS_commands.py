################################################################################
#
#  CTI_RTS_Backstop_Commands - Base Class Definition
#
################################################################################

class LTCTI_ACIS_commands(object):
    """
    The Backstop_Commands class is used when translating an ACIS
    command line found in such files as a LTCTI_CLD.txt file
    or an FOT LTCTI RTS file, into another format such as SKA.Parse
    or FOT CR*.backstop format.

         Usage: 
                 import Backstop_Commands
                 BC = Backstop_commands()
                 expanded_cmds = PCLD.ProcessCLD(CLD_file_path, CLD_START_TIME, CAP_NUMBER)

 Sample output, contents of expanded_cmds:


     NOTES:
           1) Each ACIS Ops LTCTI CLD file starts with at least a 1.025 second wait.
              Therefore the first command in the array will occur 1.025 seconds
              after the start time argument of the call to ProcessCLD.

           2) Every ACIS Ops LTCTI CLD file contains the SCS slot in which the CLD will be run.

    """

    # Constructor
    def __init__(self):

         # Initialize the recognized command list
        self.cmd_list = ['WSVIDALLDN',
                         'XTZ0000005', 
                         'XCZ0000005',
                         'RS_0000001',
                         'RH_0000001',
                         'AA00000000',
                         'WSPOW00000',
                         'WSPOW0002A',
                         'WSPOW08002',
                         'WSPOW08F3E',
                         'WSPOW08E1E',
                         'WT007AC024',
                         'WT007AE024',
                         'WT00B26014',
                         'WT00C62014',
                         'WT00C60014']
         
        self.cmd_vals = {'XTZ0000005': {'CMDS': '4',
                                        'WORDS': '4',
                                        'PACKET': 'PACKET(40)= D8000040004003A000E000400000'},
                  
                         'XCZ0000005': {'CMDS': '4',
                                        'WORDS': '4',
                                        'PACKET': 'PACKET{40)= D800004000400270010000400000'},
                  
                         'RS_0000001': {'CMDS': '3',
                                        'WORDS': '3',
                                        'PACKET': 'PACKET{40)= D80000300030042002100'},
                  
                         'RH_0000001': {'CMDS': '3',
                                        'WORDS': '3',
                                        'PACKET': 'PACKET{40)= D80000300030044002300'},
                          
                         'AA00000000': {'CMDS': '3',
                                        'WORDS': '3',
                                        'PACKET': 'PACKET{40)= D80000300030603001300'},
                          
                         'WSVIDALLDN': {'CMDS': '4',
                                        'WORDS': '5',
                                        'PACKET': 'PACKET{40)= D800005000506050020000000000'},

                         'WT007AC024': {'CMDS': '87',
                                        'WORDS': '150',
                                        'PACKET': 'PACKET{40)= D80009600962D7C00090004308AC024007A10730'},                      
                         'WT007AE024': {'CMDS': '87',
                                        'WORDS': '150',
                                        'PACKET': 'PACKET{40)= D80009600962D80000900046541E024007A65790'},
                          
                         'WT00B26014': {'CMDS': '87',
                                        'WORDS': '150',
                                        'PACKET': 'PACKET{40)= D80009600963194000900049271601400B210730'},
                         
                         'WSPOW08F3E': {'CMDS': '5',
                                        'WORDS': '7',
                                        'PACKET': 'PACKET(40)= D80000700071A9E00200000008F0001003E'},
 
                         'WSPOW0CF3F': {'CMDS': '5',
                                        'WORDS': '7',
                                        'PACKET': 'PACKET(40)= D800007000701670020000000CF0001003F'},
                          
                         'WSPOW00000': {'CMDS': '5',
                                        'WORDS': '7',
                                        'PACKET': 'PACKET(40)= D8000070007030500200000000000010000'},
                          
                         'WSPOW0002A': {'CMDS': '5',
                                        'WORDS': '7',
                                        'PACKET': 'PACKET(40)= D80000700073E800020000000000001002A'},
                          
                         'WSPOW08E1E': {'CMDS': '42',
                                        'WORDS': '73',
                                        'PACKET': 'PACKET(40)= D8000070007030500200000000000019999'},
                          
                         'WT00C62014': {'CMDS': '101',
                                        'WORDS': '102',
                                        'PACKET': 'PACKET(40)= D8000070007030500200000000000011111'},
                          
                         'WT00C60014': {'CMDS': '87',
                                        'WORDS': '150',
                                        'PACKET': 'PACKET(40)= D8000960096381600090004C357001400C6237A0'}
                                    
                                      
                        } # END Command Constants Dictionary


        # Now the LTCTI CLD files are run as an RTS so the times and
        # VCDU's in this data structure are bogus.  Meaningful times
        # will be entered in a copy of this data structure when the commands
        # are integrated into a backstop history
        # NOTE: Never write into this attribute. Always make a copy
        self.general_bs_cmds = [
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

                                ] # End of definition of self.LTCTI_bs_cmds
