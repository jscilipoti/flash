! MODULE: inputdata
! This module store some parameters and variables needed to start the flash
! calc.
module inputData
    
    use iso_fortran_env, only: int8, int32, real64
    
    implicit none

    private :: &
        & icalc, model, ipr, iout, novap, ig, &
        & name_maxlen, &
        & flash_input_filename, &
        & ntext, &
        & ant
    

    public :: &
        & nms, nmg, nint, &
        & ipareq, &
        & output, &
        & aint1

    ! Max site number per group.
    integer(kind=int32),parameter :: nms = 2
    ! Max size of the arrays which contains subgroups' data .
    integer(kind=int32),parameter :: nmg = 150
    ! Max size of the arrays which contains the interaction parameters.
    integer(kind=int32),parameter :: nint = 70

    ! The numbers on the second line of the flash-parameters-file:
    !   icalc:  [0,1,2]     (DEFAULT = 0)
    !           0-  ' **** FLASH CALCULATION ****'
    !           1-  ' **** BINODAL CURVE CALCULATION ****'
    !           2-  ' **** CALCULATION OF UNIQUAC PARAMETERS FROM UNIFAC **** '
    !
    !   model:  [0,1,2,3]   (DEFAULT = 0)
    !           0-  ' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIFAC'
    !           1-  ' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIQUAC'
    !           2-  ' MODEL USED FOR LIQUID PHASE NON-IDEALITY: A-UNIFAC'
    !           3-  ' MODEL USED FOR GROUP CONTRIBUTION MODE: GC'
    !
    !   ipr:    [0,1]       (DEFAULT = 0)
    !           1-  ' ** COMPARISON OF ACTIVITIES CALCULATED BY UNIFAC AND 
    !               UNIQUAC, RESPECTIVELY **'
    !
    !   iout:   [0,1]       (DEFAULT = 1)
    !           1-  'open 'lleasoccuzada.OUT''
    !
    !   novap:  [0,1]       (DEFAULT = 0)
    !           0-'VAPOR PHASE INCLUDED IN FLASH-CALCULATIONS'
    !
    !   ig:     [0,1]       (DEFAULT = 1)
    !           0-  'write compositions'
    !           1-  'write compositions and activities'
    !
    !   ipareq: [0,1]       (DEFAULT = 2)
    !           1-'liquid-liquid parameters table (UNIFAC)'
    !           2-'liquid-vapor parameters table (UNIFAC)'
    !           3-'infinite dilution parameters table (UNIFAC)'
    !           4-'GC-EOS parameters'
    !
    integer(kind=int8) :: &
        & icalc = 0,& 
        & model = 0,& 
        & ipr = 0,& 
        & iout = 1,& 
        & novap = 0,& 
        & ig = 1,&
        & ipareq = 2 

    integer(kind=int32) :: output
    
    ! The name of the file with the filename of the flash-parameters-file.
    !character(len=8) :: input_filename = "name.dat"
    
    ! The max character lenght of the filename with the flash parameters.
    integer(kind=int32), parameter :: name_maxlen = 36
    
    ! The filename of the flash-parameters-file.
    character(len=name_maxlen) :: flash_input_filename
    
    ! A matrix with the size of (the max site number per group x
    ! Max size of the arrays which contains subgroups' data)
    real(kind=real64), dimension(NMG,NMG) :: aint1 = 0.D0

    ! A string with the first line of the flash calc input file. It is usually
    ! the same name in the namefile
    !integer(kind=int32), dimension(name_maxlen) :: NTEXT
    !characterlen=2), dimension(name_maxlen) :: NTEXT
    character(len=name_maxlen) :: ntext
    
    real(kind=real64), dimension(10,3):: ant = 0.D0
    
    ! The original InputData module content:
    !integer,parameter::NMS = 2 !n�mero m�ximo de sitios por grupo
    !integer,parameter::NMG = 150 !(dimensi�n m�xima para vectores que guardan inf. sobre subgrupos)
    !integer,parameter::NINT = 70 !Dimensi�n para vectores que guardan inf. sobre par�metros de interacci�n
    !integer::ipareq
    !integer::output
    !real*8,dimension(NMG,NMG)::aint1

end module inputData