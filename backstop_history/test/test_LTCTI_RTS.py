import LTCTI_RTS
import LTCTI_ACIS_commands


# Create an instance of LTCTI_RTS
RTS = LTCTI_RTS.LTCTI_RTS('/home/gregg/MODELS/ACIS_HISTORY_ASSEMBLY/RTS/')

RTS.RTS_file_loc = '/home/gregg/MODELS/ACIS_HISTORY_ASSEMBLY/RTS/'

#RTS.parse_FOT_request(fot_request)
RTS.RTS_name = '1_4_CTI'
RTS.SCS_NUM = '135'
RTS.NUM_HOURS = '000:24:00.00'
RTS_start_date = '2017:251:20:45:00'

print RTS.RTS_name, RTS.SCS_NUM, RTS.NUM_HOURS, RTS_start_date

processed_commands = RTS.processRTS(RTS.RTS_name, RTS.SCS_NUM, RTS.NUM_HOURS, RTS_start_date)

ska_cmds = RTS.convert_ACIS_RTS_to_ska_parse(processed_commands)

rts_load = open(RTS.RTS_file_loc+RTS.RTS_name+'.RTS', 'r')
for eachline in rts_load:
    if (eachline[0] != '!') and (eachline[0:2] != '\n'): 
        split_line = ''.join(eachline.split())
        split_line = split_line.split(',')

        print split_line


rts_load.close()

