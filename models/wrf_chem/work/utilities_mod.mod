	  ļ4     k820309    4          12.1        ź@Q                                                                                                           
       ../../../utilities/utilities_mod.f90 UTILITIES_MOD       $       FILE_EXIST GET_UNIT OPEN_FILE CLOSE_FILE TIMESTAMP REGISTER_MODULE ERROR_HANDLER TO_UPPER NC_CHECK NEXT_FILE LOGFILEUNIT NMLFILEUNIT FIND_TEXTFILE_DIMS FILE_TO_TEXT INITIALIZE_UTILITIES FINALIZE_UTILITIES DUMP_UNIT_ATTRIBUTES FIND_NAMELIST_IN_FILE CHECK_NAMELIST_READ DO_NML_TERM SET_TASKNUM SET_OUTPUT DO_OUTPUT SET_NML_OUTPUT DO_NML_FILE E_DBG E_MSG E_WARN E_ERR DEBUG MESSAGE WARNING FATAL IS_LONGITUDE_BETWEEN GET_NEXT_FILENAME ASCII_FILE_FORMAT                      @                              
       R8 PI                      @                             
                                                                                                                                                              
                
                 -DTū!	@        3.14159265358979323846                                                                                                    0$         @                                    P                      #NCERR                      
   @                                                                                                                                             0                                             	                                          ’’’’’’’’                                                     
                                          ’’’’’’’’                                                                                                            0                                                                                                   1                                                                                                   2                                                                                                   1                                                                                                   2          @@                                                     @                                            %         @                                                          #FILE_EXIST%LEN_TRIM    #FILE_NAME                  @                                 LEN_TRIM           
  @@                                                 1 %         @                                                            %         @                                                          #OPEN_FILE%MIN    #OPEN_FILE%TRIM    #OPEN_FILE%PRESENT    #OPEN_FILE%LEN    #FNAME    #FORM    #ACTION                  @                                 MIN               @                                 TRIM               @                                 PRESENT               @                                 LEN           
  @@                                                 1           
 @@                                                 1           
 @@                                                 1 #         @                                                    #IUNIT              
   @                                         #         @                                                     #TIMESTAMP%ADJUSTL !   #TIMESTAMP%TRIM "   #STRING1 #   #STRING2 $   #STRING3 %   #POS &                 @                            !     ADJUSTL               @                            "     TRIM           
 @@                             #                    1           
 @@                             $                    1           
 @@                             %                    1           
  @@                             &                    1 #         @                                 '                  #REGISTER_MODULE%TRIM (   #SRC )   #REV *   #RDATE +                 @                            (     TRIM           
  @@                             )                    1           
  @@                             *                    1           
  @@                             +                    1 #         @                                 ,               	   #ERROR_HANDLER%TRIM -   #ERROR_HANDLER%PRESENT .   #LEVEL /   #ROUTINE 0   #TEXT 1   #SRC 2   #REV 3   #RDATE 4   #AUT 5   #TEXT2 6   #TEXT3 7                 @                            -     TRIM               @                            .     PRESENT           
   @                              /                     
  @@                             0                    1           
  @@                             1                    1           
 @@                             2                    1           
 @@                             3                    1           
 @@                             4                    1           
 @@                             5                    1           
 @@                             6                    1           
 @@                             7                    1 #         @                                 8                  #TO_UPPER%CHAR 9   #TO_UPPER%ICHAR :   #TO_UPPER%LEN ;   #STRING <                 @                            9     CHAR               @                            :     ICHAR               @                            ;     LEN           
D @@                             <                     1 #         @                                  =                  #NC_CHECK%TRIM >   #NC_CHECK%PRESENT ?   #ISTATUS @   #SUBR_NAME A   #CONTEXT B                 @                            >     TRIM               @                            ?     PRESENT           
  @@                              @                     
  @@                             A                    1           
 @@                             B                    1 $        @                                C                           #NEXT_FILE%LEN_TRIM D   #NEXT_FILE%ADJUSTL E   #NEXT_FILE%TRIM F   #NEXT_FILE%LEN G   #FNAME H   #IFILE I   H r G     5 O p                                  @                            D     LEN_TRIM               @                            E     ADJUSTL               @                            F     TRIM               @                            G     LEN           
  @@                             H                    1           
   @                              I           #         @                                 J                  #FIND_TEXTFILE_DIMS%LEN_TRIM K   #FIND_TEXTFILE_DIMS%TRIM L   #FIND_TEXTFILE_DIMS%PRESENT M   #FNAME N   #NLINES O   #LINELEN P                 @                            K     LEN_TRIM               @                            L     TRIM               @                            M     PRESENT           
  @@                             N                    1           D  @                              O                      F @@                              P            #         @                                  Q                  #FILE_TO_TEXT%SIZE R   #FILE_TO_TEXT%MIN S   #FILE_TO_TEXT%TRIM T   #FILE_TO_TEXT%LEN U   #FNAME V   #TEXTBLOCK W                 @                            R     SIZE               @                            S     MIN               @                            T     TRIM               @                            U     LEN           
  @@                             V                    1 ,          D@@                             W                                   &                                           1 #         @                                 X                  #INITIALIZE_UTILITIES%ADJUSTL Y   #INITIALIZE_UTILITIES%TRIM Z   #INITIALIZE_UTILITIES%PRESENT [   #PROGNAME \   #ALTERNATENAME ]   #OUTPUT_FLAG ^                 @                            Y     ADJUSTL               @                            Z     TRIM               @                            [     PRESENT           
 @@                             \                    1           
 @@                             ]                    1           
 @@                              ^           #         @                                 _                  #FINALIZE_UTILITIES%TRIM `   #FINALIZE_UTILITIES%PRESENT a   #PROGNAME b                 @                            `     TRIM               @                            a     PRESENT           
 @@                             b                    1 #         @                                  c                  #DUMP_UNIT_ATTRIBUTES%ADJUSTL d   #DUMP_UNIT_ATTRIBUTES%TRIM e   #IUNIT f                 @                            d     ADJUSTL               @                            e     TRIM           
   @                              f           #         @                                 g                  #FIND_NAMELIST_IN_FILE%ADJUSTL h   #FIND_NAMELIST_IN_FILE%TRIM i   #FIND_NAMELIST_IN_FILE%PRESENT j   #NAMELIST_FILE_NAME k   #NML_NAME l   #IUNIT m   #WRITE_TO_LOGFILE_IN n                 @                            h     ADJUSTL               @                            i     TRIM               @                            j     PRESENT           
  @@                             k                    1           
  @@                             l                    1           D  @                              m                      
 @@                              n           #         @                                 o                  #CHECK_NAMELIST_READ%INDEX p   #CHECK_NAMELIST_READ%TRIM q   #CHECK_NAMELIST_READ%PRESENT r   #CHECK_NAMELIST_READ%LEN s   #IUNIT t   #IOSTAT_IN u   #NML_NAME v   #WRITE_TO_LOGFILE_IN w                 @                            p     INDEX               @                            q     TRIM               @                            r     PRESENT               @                            s     LEN           
  @@                              t                     
   @                              u                     
  @@                             v                    1           
 @@                              w           %         @                                x                            #         @                                  y                   #TASKNUM z             
   @                              z           #         @                                  {                   #DOFLAG |             
   @                              |           %         @                                }                            #         @                                 ~                  #SET_NML_OUTPUT%TRIM    #NMLSTRING                  @                                 TRIM           
  @@                                                 1 %         @                                                            %         @                                                           #IS_LONGITUDE_BETWEEN%MODULO    #IS_LONGITUDE_BETWEEN%PRESENT    #LON    #MINLON    #MAXLON    #DORADIANS                  @                                 MODULO               @                                 PRESENT           
  @@                                  
                
  @@                                  
                
  @@                                  
                
 @@                                         $         @                                                          #GET_NEXT_FILENAME%LEN_TRIM    #GET_NEXT_FILENAME%ADJUSTL    #GET_NEXT_FILENAME%LEN    #LISTNAME    #INDEX                          @                                 LEN_TRIM               @                                 ADJUSTL               @                                 LEN           
  @@                                                 1           
   @                                         %         @                                                           #ASCII_FILE_FORMAT%ADJUSTL    #ASCII_FILE_FORMAT%TRIM    #ASCII_FILE_FORMAT%PRESENT    #ASCII_FILE_FORMAT%LEN    #FFORM                  @                                 ADJUSTL               @                                 TRIM               @                                 PRESENT               @                                 LEN           
 @@                                                 1        ;      fn#fn #   Ū   Ņ  b   uapp(UTILITIES_MOD    ­  F   J  TYPES_MOD    ó  @   J  NETCDF    3  p       R8+TYPES_MOD    £         PI+TYPES_MOD "   )  q       NF90_NOERR+NETCDF %     c       NF90_STRERROR+NETCDF +   ż  @   e   NF90_STRERROR%NCERR+NETCDF    =  q       MESSAGE    ®  p       DEBUG      p       E_DBG      q       E_MSG    ’  q       E_WARN    p  q       E_ERR    į  q       WARNING    R  q       FATAL    Ć  @       LOGFILEUNIT    	  @       NMLFILEUNIT    C	  x       FILE_EXIST $   »	  A      FILE_EXIST%LEN_TRIM %   ü	  L   a   FILE_EXIST%FILE_NAME    H
  P       GET_UNIT    
  Ā       OPEN_FILE    Z  <      OPEN_FILE%MIN      =      OPEN_FILE%TRIM "   Ó  @      OPEN_FILE%PRESENT      <      OPEN_FILE%LEN     O  L   a   OPEN_FILE%FNAME      L   a   OPEN_FILE%FORM !   ē  L   a   OPEN_FILE%ACTION    3  S       CLOSE_FILE !     @   a   CLOSE_FILE%IUNIT    Ę  £       TIMESTAMP "   i  @      TIMESTAMP%ADJUSTL    ©  =      TIMESTAMP%TRIM "   ę  L   a   TIMESTAMP%STRING1 "   2  L   a   TIMESTAMP%STRING2 "   ~  L   a   TIMESTAMP%STRING3    Ź  L   a   TIMESTAMP%POS              REGISTER_MODULE %     =      REGISTER_MODULE%TRIM $   Ņ  L   a   REGISTER_MODULE%SRC $     L   a   REGISTER_MODULE%REV &   j  L   a   REGISTER_MODULE%RDATE    ¶  Ł       ERROR_HANDLER #     =      ERROR_HANDLER%TRIM &   Ģ  @      ERROR_HANDLER%PRESENT $     @   a   ERROR_HANDLER%LEVEL &   L  L   a   ERROR_HANDLER%ROUTINE #     L   a   ERROR_HANDLER%TEXT "   ä  L   a   ERROR_HANDLER%SRC "   0  L   a   ERROR_HANDLER%REV $   |  L   a   ERROR_HANDLER%RDATE "   Č  L   a   ERROR_HANDLER%AUT $     L   a   ERROR_HANDLER%TEXT2 $   `  L   a   ERROR_HANDLER%TEXT3    ¬         TO_UPPER    9  =      TO_UPPER%CHAR    v  >      TO_UPPER%ICHAR    “  <      TO_UPPER%LEN     š  L   a   TO_UPPER%STRING    <         NC_CHECK    Ö  =      NC_CHECK%TRIM !     @      NC_CHECK%PRESENT !   S  @   a   NC_CHECK%ISTATUS #     L   a   NC_CHECK%SUBR_NAME !   ß  L   a   NC_CHECK%CONTEXT    +  š       NEXT_FILE #     A      NEXT_FILE%LEN_TRIM "   \  @      NEXT_FILE%ADJUSTL      =      NEXT_FILE%TRIM    Ł  <      NEXT_FILE%LEN       L   a   NEXT_FILE%FNAME     a  @   a   NEXT_FILE%IFILE #   ”  Ź       FIND_TEXTFILE_DIMS ,   k  A      FIND_TEXTFILE_DIMS%LEN_TRIM (   ¬  =      FIND_TEXTFILE_DIMS%TRIM +   é  @      FIND_TEXTFILE_DIMS%PRESENT )   )  L   a   FIND_TEXTFILE_DIMS%FNAME *   u  @   a   FIND_TEXTFILE_DIMS%NLINES +   µ  @   a   FIND_TEXTFILE_DIMS%LINELEN    õ  ¼       FILE_TO_TEXT "   ±  =      FILE_TO_TEXT%SIZE !   ī  <      FILE_TO_TEXT%MIN "   *  =      FILE_TO_TEXT%TRIM !   g  <      FILE_TO_TEXT%LEN #   £  L   a   FILE_TO_TEXT%FNAME '   ļ     a   FILE_TO_TEXT%TEXTBLOCK %      Ż       INITIALIZE_UTILITIES -   \!  @      INITIALIZE_UTILITIES%ADJUSTL *   !  =      INITIALIZE_UTILITIES%TRIM -   Ł!  @      INITIALIZE_UTILITIES%PRESENT .   "  L   a   INITIALIZE_UTILITIES%PROGNAME 3   e"  L   a   INITIALIZE_UTILITIES%ALTERNATENAME 1   ±"  @   a   INITIALIZE_UTILITIES%OUTPUT_FLAG #   ń"         FINALIZE_UTILITIES (   #  =      FINALIZE_UTILITIES%TRIM +   Į#  @      FINALIZE_UTILITIES%PRESENT ,   $  L   a   FINALIZE_UTILITIES%PROGNAME %   M$         DUMP_UNIT_ATTRIBUTES -   į$  @      DUMP_UNIT_ATTRIBUTES%ADJUSTL *   !%  =      DUMP_UNIT_ATTRIBUTES%TRIM +   ^%  @   a   DUMP_UNIT_ATTRIBUTES%IUNIT &   %  ų       FIND_NAMELIST_IN_FILE .   &  @      FIND_NAMELIST_IN_FILE%ADJUSTL +   Ö&  =      FIND_NAMELIST_IN_FILE%TRIM .   '  @      FIND_NAMELIST_IN_FILE%PRESENT 9   S'  L   a   FIND_NAMELIST_IN_FILE%NAMELIST_FILE_NAME /   '  L   a   FIND_NAMELIST_IN_FILE%NML_NAME ,   ė'  @   a   FIND_NAMELIST_IN_FILE%IUNIT :   +(  @   a   FIND_NAMELIST_IN_FILE%WRITE_TO_LOGFILE_IN $   k(        CHECK_NAMELIST_READ *   o)  >      CHECK_NAMELIST_READ%INDEX )   ­)  =      CHECK_NAMELIST_READ%TRIM ,   ź)  @      CHECK_NAMELIST_READ%PRESENT (   **  <      CHECK_NAMELIST_READ%LEN *   f*  @   a   CHECK_NAMELIST_READ%IUNIT .   ¦*  @   a   CHECK_NAMELIST_READ%IOSTAT_IN -   ę*  L   a   CHECK_NAMELIST_READ%NML_NAME 8   2+  @   a   CHECK_NAMELIST_READ%WRITE_TO_LOGFILE_IN    r+  P       DO_NML_TERM    Ā+  U       SET_TASKNUM $   ,  @   a   SET_TASKNUM%TASKNUM    W,  T       SET_OUTPUT "   «,  @   a   SET_OUTPUT%DOFLAG    ė,  P       DO_OUTPUT    ;-  p       SET_NML_OUTPUT $   «-  =      SET_NML_OUTPUT%TRIM )   č-  L   a   SET_NML_OUTPUT%NMLSTRING    4.  P       DO_NML_FILE %   .  Ć       IS_LONGITUDE_BETWEEN ,   G/  ?      IS_LONGITUDE_BETWEEN%MODULO -   /  @      IS_LONGITUDE_BETWEEN%PRESENT )   Ę/  @   a   IS_LONGITUDE_BETWEEN%LON ,   0  @   a   IS_LONGITUDE_BETWEEN%MINLON ,   F0  @   a   IS_LONGITUDE_BETWEEN%MAXLON /   0  @   a   IS_LONGITUDE_BETWEEN%DORADIANS "   Ę0  Ė       GET_NEXT_FILENAME +   1  A      GET_NEXT_FILENAME%LEN_TRIM *   Ņ1  @      GET_NEXT_FILENAME%ADJUSTL &   2  <      GET_NEXT_FILENAME%LEN +   N2  L   a   GET_NEXT_FILENAME%LISTNAME (   2  @   a   GET_NEXT_FILENAME%INDEX "   Ś2  Š       ASCII_FILE_FORMAT *   Ŗ3  @      ASCII_FILE_FORMAT%ADJUSTL '   ź3  =      ASCII_FILE_FORMAT%TRIM *   '4  @      ASCII_FILE_FORMAT%PRESENT &   g4  <      ASCII_FILE_FORMAT%LEN (   £4  L   a   ASCII_FILE_FORMAT%FFORM 