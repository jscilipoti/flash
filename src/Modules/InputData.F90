module InputData
    
    integer,parameter::NMS = 2 !n�mero m�ximo de sitios por grupo
    integer,parameter::NMG = 150 !(dimensi�n m�xima para vectores que guardan inf. sobre subgrupos)
    integer,parameter::NINT = 70 !Dimensi�n para vectores que guardan inf. sobre par�metros de interacci�n
    integer::ipareq
    integer::output
    real*8,dimension(NMG,NMG)::aint1

endmodule InputData