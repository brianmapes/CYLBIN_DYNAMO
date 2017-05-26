
IDL_CYLBIN_process_TOGA.pro:

It is the main IDL program where CYLBIN processing is done. Subscript "TOGA" is added because, there are 
two version of the main program that is being used to process both TOGA radar data from the revelle ship platform 
and SMARTR data from the gan island. The difference between the two versions is only some of the lines are commented 
in the IDL main program that is being used for the SMARTR data as the variable "SDZhsit” was not available for the SMARTR data.

IDL_CYLBIN_process_SMARTR.pro

It is the main IDL program where CYLBIN processing is done for the SMARTR data during the DYNAMO. The difference 
between the two versions (-TOGA and -SMARTR) of the program are explained above.

 
IDL_CYLBIN_params.pro

It is the IDL function that is called from inside the main program, where the directory path for the input/output files, and the parameter 
used in the main IDL program are defined.


IDL_CYLBIN_proc_subroutines.pro

It is the IDL function that is called from inside the main program, where different 
tasks are performed as different modules.


