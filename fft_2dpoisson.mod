  �/  ^   k820309    p          19.1        �l_                                                                                                          
       fft_2dpoisson.f90 FFT_2DPOISSON                                                                                                     0                                                                                                   1                                                                                                   2                                                                                                   3                                                                                                   4                                                                                                   5                                                                                                   6                                                                                                   7                                             	                                                      8                                             
                                       	               9                                                                                    
               10                                                                                       ��������                                                                                                           1                                                                                                    0                                                                                                   1                                                                                                   2                                                                                                   4                                                                                                   8                                                                                                   16                                                                                                    32                                                                                    @               64                                                                                                    2097152                                                                                    �               128                                                                                                   256                                                                                                   512                                                                                                   1024                                                                                                   2048                                                                                                   4096                                                                                                    8192                                                                                     @              16384                                                                                     �              32768                                                                                                    65536                                             !                                                      131072                                             "                                                      262144                                             #                                                      524288                                             $                                                      1048576                                             %     
                   
                  -DT�!	@                                                    &     ACOS                                              '                                                      1                                             (                                                       32                                             )                                                      1                                             *                                                       32                                             +     
                   
                  -DT�!@                                                     ,     
                   
                  -DT�!@                                                     -     
                 
                       �?        1.0D0           @@                               .                       @                                /     
                  @                                0     
                  @                                1     
                  @                                2     
                  @                                3     
                  @                                4     
                  @                                5     
                  @                                6     
                  @                                7     
                  @    �                            8     �             
      p           & p         p !         & p         p !           p "         p "                                    @    �                            9     �             
      p           & p         p !         & p         p !           p "         p "                         #         @                                   :                     #         @                                   ;                    #VX_P <   #VY_P =             D     �                            <     A             
     p           & p         p           & p         p             p !         p !                                   D     �                            =     A             
     p           & p         p           & p         p             p !         p !                         #         @                                   >                    #PHI_TH ?             D     �                            ?                   
     p !         & p        p           & p        p             p           p                           #         @                                  @                    #VX_P A   #VY_P B   #DIVU C             
      �                            A     A             
    p           & p         p           & p         p             p !         p !                                   
      �                            B     A             
    p           & p         p           & p         p             p !         p !                                   D     �                            C                   
     p !         & p        p           & p        p             p           p                           #         @                                  D                    #RHS E   #RHSF F             
@ @   �                            E                    
 	   p          & p        p             p                                     D @   �                           F                    
    p           & p         p            p                          #         @                                  G                    #PF H   #P I             
@ @   �                           H                       p           & p         p            p                                    D @   �                            I                    
     p          & p        p             p                           #         @                                  J                    #DIVU K   #DIVUF L             
      �                            K                   
    p !         & p        p           & p        p             p           p                                     D     �                           L                        p          & p         p          & p        p             p          p                           #         @                                  M                    #PHIF N   #PHI O             
      �                           N     B                 p           & p         p          & p         p !           p          p "                                   D     �                            O     �             
     p           & p         p !         & p         p !           p "         p "                         #         @                                  P                    #PHIF Q             D     �                           Q     B                  p           & p         p          & p         p !           p          p "                         #         @                                  R                    #PHI S             D     �                            S     �             
     p           & p         p !         & p         p !           p "         p "                         #         @                                  T                    #DIVUF U   #PHIF V             
      �                           U                       p          & p         p          & p        p             p          p                                     D @   �                           V     B                  p           & p         p          & p         p !           p          p "                         #         @                                   W                    #VX_P X   #VY_P Y   #PHI Z             
  @   �                            X     A             
    p           & p         p           & p         p             p !         p !                                   
  @   �                            Y     A             
    p           & p         p           & p         p             p !         p !                                   D @   �                            Z     �             
     p           & p         p !         & p         p !           p "         p "                         #         @                                   [                    #PHI \   #PHI_TH ]             
      �                            \     �             
     p           & p         p !         & p         p !           p "         p "                                   
      �                            ]                   
 !   p !         & p        p           & p        p             p           p                              �   (      fn#fn    �   q       FFTW_R2HC    9  q       FFTW_HC2R    �  q       FFTW_DHT      q       FFTW_REDFT00    �  q       FFTW_REDFT01    �  q       FFTW_REDFT10    n  q       FFTW_REDFT11    �  q       FFTW_RODFT00    P  q       FFTW_RODFT01    �  q       FFTW_RODFT10    2  r       FFTW_RODFT11    �  p       FFTW_FORWARD      q       FFTW_BACKWARD    �  q       FFTW_MEASURE #   �  q       FFTW_DESTROY_INPUT    g  q       FFTW_UNALIGNED %   �  q       FFTW_CONSERVE_MEMORY     I  q       FFTW_EXHAUSTIVE $   �  r       FFTW_PRESERVE_INPUT    ,	  r       FFTW_PATIENT    �	  r       FFTW_ESTIMATE !   
  w       FFTW_WISDOM_ONLY &   �
  s       FFTW_ESTIMATE_PATIENT #   �
  s       FFTW_BELIEVE_PCOST !   m  s       FFTW_NO_DFT_R2HC $   �  t       FFTW_NO_NONTHREADED "   T  t       FFTW_NO_BUFFERING $   �  t       FFTW_NO_INDIRECT_OP )   <  t       FFTW_ALLOW_LARGE_GENERIC $   �  u       FFTW_NO_RANK_SPLITS %   %  u       FFTW_NO_VRANK_SPLITS !   �  u       FFTW_NO_VRECURSE      v       FFTW_NO_SIMD    �  v       FFTW_NO_SLOW ,   �  v       FFTW_NO_FIXED_RADIX_LARGE_N #   q  w       FFTW_ALLOW_PRUNING    �  p       PI    X  =       ACOS    �  q       NXMIN      r       NXMAX    x  q       NYMIN    �  r       NYMAX    [  p       XMAX    �  p       YMAX    ;  u       DT    �  @       NG    �  @       DX    0  @       DY    p  @       DDT    �  @       DDX    �  @       DDY    0  @       DDX2    p  @       DDY2    �  @       XMIN    �  @       YMIN    0  �       X      �       Y    �  H       SET_GRID       \       SET_VELOCITY "   |  �   a   SET_VELOCITY%VX_P "   P  �   a   SET_VELOCITY%VY_P    $  T       CAL_THEORY "   x  �   a   CAL_THEORY%PHI_TH     L  f       CAL_POISSON_RHS %   �  �   a   CAL_POISSON_RHS%VX_P %   �  �   a   CAL_POISSON_RHS%VY_P %   Z  �   a   CAL_POISSON_RHS%DIVU    .  [       FFT_1D_EXE    �  �   a   FFT_1D_EXE%RHS     -   �   a   FFT_1D_EXE%RHSF    �   W       IFFT_1D_EXE    (!  �   a   IFFT_1D_EXE%PF    �!  �   a   IFFT_1D_EXE%P    p"  ]       FFT_RHS    �"  �   a   FFT_RHS%DIVU    �#  �   a   FFT_RHS%DIVUF    u$  [       IFFT_PHIF    �$  �   a   IFFT_PHIF%PHIF    �%  �   a   IFFT_PHIF%PHI    x&  R       SET_BC_ITR     �&  �   a   SET_BC_ITR%PHIF    �'  Q       SET_BC    �'  �   a   SET_BC%PHI    �(  ]       ITR_POISSON "    )  �   a   ITR_POISSON%DIVUF !   �)  �   a   ITR_POISSON%PHIF    �*  e       CAL_POISSON !   -+  �   a   CAL_POISSON%VX_P !   ,  �   a   CAL_POISSON%VY_P     �,  �   a   CAL_POISSON%PHI    �-  ]       CAL_ERROR    .  �   a   CAL_ERROR%PHI !   �.  �   a   CAL_ERROR%PHI_TH 