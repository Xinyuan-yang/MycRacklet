cRacklet: a spectral boundary integral method library for interfacial rupture simulation
===========================================================================================

cRacklet is a C++ boundary integral library based on a spectral formulation of the dynamic wave
equations in two semi-infinite linearly elastic solids. Its implementation is specially tailored
for the modeling of dynamic crack/rupture propagation along a planar interface bonding the two
half-spaces. The main benefit of this spectral method is a numerical discretization limited to
the mid-plane, thereby providing a very fine description of the dynamic rupture processes,
unattainable with finite-element or finite-difference schemes.
For more details about the method and its applications, refer to "related publications" section.

## Dependencies

The following dependencies are required to compile cRacklet:

- a C++ compiler
- [cMake](https://cmake.org/)
- [FFTW3](http://www.fftw.org/)
- [GSL](https://www.gnu.org/software/gsl/)

Optional dependencies are:

- a Fortran compiler to generate new kernels
- [python 3+](https://www.python.org/), for python binding.
- [pybind11](https://github.com/pybind/pybind11), for python binding, automatically installed if not found on the system.
- FFTW3 with OpenMP support, for multi-threaded parallel computing
- [pytest](https://docs.pytest.org/en/latest/) (for tests)
- [Doxygen](http://doxygen.nl/) (for documentation)
- [Sphinx](https://www.sphinx-doc.org/en/stable/) (for documentation)

## Compile and build

You can compile cRacklet using cMake:

    mkdir build
    cd build
    ccmake .. or cmake ..
    make

## Tests

You need to activate the test options with cMake (`CRACKLET_TESTS`). Go inside `build/` and then run the test with:

    make test

## Documentation

Online documentation is available at https://cracklet.gitlab.io/cracklet/ .

You can generate the documentation locally using sphinx and doxygen. You first have to activate the option `CRACKLET_DOCUMENTATION` with cmake. Then go inside `build/` and run

    make dev_doc

The documentation will be generated in html format and store inside build/doc/html

## Usage - Examples

Exemples of simulations are available in the folder `examples`.

You can create your own CMake project and link it to cRacklet by adding

    find_package(cRacklet REQUIRED)

to your `CMakeLists`. After doing so, you can create an executable with the command

    add_cracklet_simulation(executable source_file
    NU_TOP nu_t
    NU_BOT nu_b)

with `NU_TOP` and `NU_BOT` optional arguments that, if provided, will copy the corresponding kernels to your simulation folder. Note that the convolution kernels describing the elastodynamic response
of the surrounding solids (as function of Poisson's ratio) shall be pre-computed using the Fortran routine
`inverse_serial.f` provided in the folder `pre-computed-kernels/`. The kernels for $`\nu = 0.33`$ and $`\nu = 0.35`$ are available by default.

Basically, the cRacklet engine is composed of seven major types of object:

### SpectralModel:   
This class is the core of cRacklet librairy and contains the methods processing the different steps
required to solve the elastodynamic response of the two semi-infinite half space.

### InterfaceLaw:
This abstract class contains the law describing interface conditions (e.g. cohesive law, frictional interface).
New interface laws can be simply added by creating classes inheriting from InterfaceLaw framework. 

### ContactLaw:
This abstract class contains law describing contact law in case overlapping is prevented at the interface.
New contact laws can be added by creating classes inheriting from ContactLaw.

### Interfacer:
This object helps the user to create and define initial interface properties.
Interfacer is templated by the chosen type of InterfaceLaw.

### DataRegister:
This class is the static centralized memory register which contains most of the simulation data

### SimulationDriver:
High-level class providing a UI to drive simulations in different standard situations.

### DataDumper: 
This class is used to generate various types of output file during simulation.

## PYTHON INTERFACE

In order to use the python interface build for cRacklet, you need pybind11 and python3.

During configuration, activate `CRACKLET_PYTHON_INTERFACE`. The make command will create a python library in the `build/python/` folder. Please add this path to your python path.

You can activate the python example with the option `CRACKLET_EXAMPLES_PYTHON`

An example of the use of the python interface is provided in `build/examples/python/001-modeI/modeI.py`

## Tutorials

Tutorials with the python interface are also available on Binder with a pre-installed version of cRacklet.

[Supershear transition along heterogeneous interface](https://mybinder.org/v2/gl/cracklet%2Ftutorials/72d3a17295b091c74cd57fc5647e9e91dcfc5e51?filepath=supershear%2Fsupershear.ipynb)

[Frictional crack nucleation and propagation for rate and state friction](https://mybinder.org/v2/gl/cracklet%2Ftutorials/2cbe12707784f451c14ffe401e95ab5d6b415b8b?filepath=rate-and-state%2Frate_and_state.ipynb)

The source for the notebook tutorials is available [here](https://gitlab.com/cracklet/tutorials)

## Contributing

Contributions to cRacklet are welcome! Please follow the guidelines below.

### Report an issue

If you have an account on [gitlab](https://gitlab.com), you can [submit an
issue](https://gitlab.com/cracklet/cracklet/-/issues/new). The full list of issues
is available [here](https://gitlab.com/cracklet/cracklet/-/issues).

### Submit a patch / merge-request

Follow [this
guide](https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html#new-merge-request-from-a-fork)
to create a merge request on GitLab. Please target the repository's `develop`
branch.


## AUTHORS

- Fabian Barras <fabian.barras@epfl.ch>
- Thibault Roch <thibault.roch@epfl.ch>
- Damien Spielmann
- David Kammer
- Guillaume Anciaux
- Nicolas Richart
- Philippe H Geubelle
- Jean-François Molinari

## License

cRacklet is distributed under the terms of the [GNU General Public License
v3.0](https://www.gnu.org/licenses/gpl.html).

## RELATED PUBLICATIONS

The following publications have been made possible with cRacklet:

- [Barras, F., Kammer, D. S., Geubelle, P. H., Molinari, J.-F. (2014) A study of frictional contact in dynamic fracture along bimaterial interfaces. International Journal of Fracture 189(2), 149–162](https://doi.org/10.1007/s10704-014-9967-z)

- [Barras, F., Geubelle, P. H., Molinari, J.-F. (2017) Interplay between Process Zone and Material Heterogeneities for Dynamic Cracks. Physical Review Letters 119(14)](https://doi.org/10.1103/PhysRevLett.119.144101)

- [Barras, F., Carpaij, R., Geubelle, P. H., & Molinari, J.-F. (2018). Supershear bursts in the propagation of a tensile crack in linear elastic material. Physical Review E, 98(6), 063002](https://doi.org/10.1103/PhysRevE.98.063002)

- [Brener, E. A., Aldam, M., Barras, F., Molinari, J.-F., & Bouchbinder, E. (2018). Unstable Slip Pulses and Earthquake Nucleation as a Nonequilibrium First-Order Phase Transition. Physical Review Letters, 121(23), 234302](https://doi.org/10.1103/PhysRevLett.121.234302)

- [Barras, F., Aldam, M., Roch, T., Brener, E. A., Bouchbinder, E., & Molinari, J.-F. (2019). Emergence of Cracklike Behavior of Frictional Rupture: The Origin of Stress Drops. Physical Review X, 9(4), 041043](https://doi.org/10.1103/PhysRevX.9.041043)

- [Barras, F., Aldam, M., Roch, T., Brener, E. A., Bouchbinder, E., & Molinari, J.-F. (2020). The emergence of crack-like behavior of frictional rupture: Edge singularity and energy balance. Earth and Planetary Science Letters, 531, 115978](https://doi.org/10.1016/j.epsl.2019.115978)

- [Fekak, F., Barras, F., Dubois, A., Spielmann, D., Bonamy, D., Geubelle, P. H., & Molinari J. F. (2020). Crack front waves: A 3D dynamic response to a local perturbation of tensile and shear cracks. Journal of the Mechanics and Physics of Solids, 135, 103806](https://doi.org/10.1016/j.jmps.2019.103806)

- [Rezakhani, R., Barras, F., Brun, M., & Molinari, J.-F. (2020). Finite element modeling of dynamic frictional rupture with rate and state friction. Journal of the Mechanics and Physics of Solids, 141, 103967](https://doi.org/10.1016/j.jmps.2020.103967)

- [Brener, E. A., & Bouchbinder, E. (2021). Unconventional singularities and energy balance in frictional rupture. Nature Communications, 12(1), 2585](https://doi.org/10.1038/s41467-021-22806-9)

- [Lebihain, M., Roch, T., Violay, M., & Molinari, J.-F. (2021). Instability regimes in the onset of motion along disordered frictional surfaces. arXiv:2102.10870 [cond-Mat]](http://arxiv.org/abs/2102.10870)

- [Roch, T., Brener, E. A., Molinari, J.-F., & Bouchbinder, E. (2021). Velocity-driven frictional sliding: Coarsening and steady-state pulse trains. arXiv:2104.13110 [cond-Mat,Physics:nlin, Physics:physics]](http://arxiv.org/abs/2104.13110)

                                                                                                  
                        ´.-,,,...´´                                                                     
                       .,,,,,,,-----,,,....´´                                                           
                     ´.,,,,,,,,,,,,,,,,,,,,,,,,,....´´´                                                 
                    .,.....,.,,,,,,,,,,,,,,,,,,,,,,,,,,,,.....´´´                                       
                   ´}}[(l>+:°¹,,..,,,,,,,,,,,-,,,,,,,,,,,,,,,,,,¹,......´´                              
                   ´zhaaaaabttahh})i+~¹-,,.,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,....´´´                    
                   ´{bbahaabbbttttttmmmtau}?<+;~¹-,..,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,....´´           
                    }abaaaaaabbbtttttttttbbaaahhhau{(!<+:°¹,,.,,,,,,,,,,,,,,,,,,,,,,,,,,,---¹)          
                    ]aaaahaaaabbbttmmmtttbbahhhhaaaaaaaaaabahz]?i+:¹,,......,,,,,,,,,,,,---°OY          
                    )aaaaaa>¹++abbbtttttttbbaaaaaaabaaaaabbabbttmmYVVYbz[?<;°-,...,,,,,,,-°QBY          
                    !haaa{-....-)bbttttttbbaaaaaaabbbbbbbbbbbbttmtmYYYYVVOOOOQOYaz[?i+:°,~MWWY          
                    ¹;>l)-.......~utbttttbbbaaaaaaaaaabbbbbbbbbbtttmmYYVVVOOOOOQQQMRXXXBXRBWXY          
                         ....´.´´,.<bbttttbbbbabbaabbbbbbbbtttttbttttmmYYVVVOOOOOQZZRXWBBXBQXY          
                        ´;,´.´´´´...¹:+![hbbbbbbbbbbbbbbbtttttttttttttmYYVVVOOOOQQZZRWWBBXBmRt          
                        .++¹....´´´...´   ´.,°;>![ubbbtttttmtmmtttttttttmYYYVVOOOOZZMXWBBWWMMa          
                        .+<>~..´´´´....´           ´.-°;i){bmmmmmmmmmYmmmmYYVYVVOOQQZMXWBXXZZi          
                        ,+i<<;,´´´.´,°,..´                  .,°;>?}hYYmtmmmmYYVVVOOOQZMXBXRZz           
                        ¹<lll<+°.....¹¹´´..                         ´.¹~+!]utYYVVVOOQZZMXRMM]           
                        ¹>ii><<i;,....´.¹~...                               ´.-~+?{bVQQZRXQV[           
                        ,>i>~;+<i>°.´....°,...                                      +QQZZRYV[           
                        .><;:~;+ii<;,...°-..´.´.                                    +VOOQQYV}           
                        ´+>::~;++>ll>°..,.....,,..                                  ;VOOQQYV}           
                         :i+;~;+i<l!!i;,.....,.....´                                ;VOOOOYV{           
                         °i>;;++>lll!ll>~....,,.,,.,.                               ;OOOOOYVz           
                         ,i>>++illl!!+ii<+,...-,-..,.´´                             :OOOOOYYz           
                         .<>i<iiilli+i<ili>~,...¹.,..,,.                            ~OOQOQVVu           
                          +<llillli<;~+illli;,..,...,¹,.,´                          ~OQQQQVYh           
                        .i]ll!il!!!l+:++>i!li>°....,°,.¹,,,´                        °QQQQQVYh           
                       ,)tM?!lii!!l!i++<<l!llli:,....,.,,,°,.                       °QQQQQOYa           
                       mVZBu)(?l<i>ii>+iiilll<ll+¹....,,.,,¹-.´                     ¹QQQZQOVa           
                       RQQWR}ut[(?!i+~:;>lllii>ll<;,..,-.,.,;,,.                    ¹QQQZZOVt           
                       iZZBmWXOz]?ll>;;~il<<>+>!li>¹..,,,,,°,,,.´                   -ZZQZZOVt           
                        mZMR##WMb[)ll>;+:;+><>+ll!li:,..,,¹:,,,,,´                  ,ZZZZZQVt           
                         ;QMRR$#$RV}(i!l<<;>+i>>:>+<ll>°...,..,,,,,.                ,ZZZZZQVm           
                          aZMRRB##XQu())!!l<;:;~~:+;;lli;-....,,,,,,,´              ,MMZZZQYY           
                          ¹MMRWBX##BMb}])?!+;;::~~;;:+ill>:,..,,,.¹°,-´             .MMMMMZmm´          
                           (MMRXBW$#$RO{[))!+;;;;;::;+;+ill+¹,..,,-¹,°,´            .MMMMMZmm´          
                           .tMMMXBBB@#XZu[?ll;;;;:~:;+>><+ll<:,.,,,--¹¹.            .MMMMMZmm´          
                           °MRMMRXXBB##BMt[)l<+;+:~::+<li+>!?!+°,,,,,--.            .RMMMMMtt´          
                           ,R#MMRXRXBBBXXRO{(l!!i++;;>ilil+<!??l;-,,---´            .MRRRRMbt´          
                          ,[BBMRRRZMXW#XB#XZh[l)?!!!!?)?li<i!))?l+°,,¹.             ´MRRRRRbh.          
                         ;QQ;+i[VB#$B$$@BB#BMm{])!??))))))!??))?il>::¹              ´MRRRRRbb.          
                         [Omtmtat@YM#W$BBBB##RQu[()]((((((??!!?llli>°               ´MRRRRXtu~          
                        +ttbbbttO#??RBWRXRXW$#WZb{}[[]]((((!ll!ill+¹                ´MRXXRXYmi          
                      ´!tbbbbbbbb#YVW#WRRMRXBB#$MYu{[((()()?!l!??+.                  MXXXXXRh?          
                     .(btbaaahaaam$XRRbOMMXXBBW$#RQh}}](()????li°´                   MXXXXXRX<          
                    .}bthhahhahhhhattuhumRRXBBWXW#WYz{}[(((((l;.                     MXXXXWRt,          
                   -uhhhhahuhuuhuuhuuuuuzhOXBWWRXX$bzz}}}[(([}bu?>~-.                ZWXXXWV{,          
                  ~hhhuuzzuzzuzhzu{{{{zzzz{mRXWWWWRMRZQQMYuuuhhhhhbhhz(i;¹,´         ZWWWXWmb,          
                 +huu{z{{z}{z{}{}{{{z{{u}uzhhZWB$$$WR$$#Mzhuhuhhuuhuhhhzuuzz(l;~,.   ZWWWWWYa°.´        
               ´lzu}uzz}}{}{}{{[}{[}[}z}z}u{zzVOM$#$BMRZa}{{uzzuzzzzzz{zz{{{{{uzuu{)iRWWWWBVOYBO,       
              .(zzz{}{{z}}}[]{[}[}}}{{}}[}}[{{z!)}zaZYbt}z{}{}{}}{{[}}}}}}{}}{}}}{uzzXBBWWBR$@#W+       
             .}{{}[{}{]][][[{][(]}]][[)]}]{}{[h()([}]}[}}{[}][[}}]}{}{[}[}[}[{{[z{z}{RWWWWBXB@RX;       
            ¹}u{{}}[[][[}](([)[((](([[(](]((()h]!]((}][]](]][[}}}[[]}[[[[[][]][[[[[{zRBBWWBX#@BX¹       
           ~uz{[[[[]]((()(](]](?((](}]([?[[?)!(];(])]][[(](](([][[[([[]](([())[[]([[}ZBBBWWZV#$X.       
          ;zz[](}(??]()(])()]))[][]}[[z][(]}[]<:<[](!!)()(l?]()((()ll!)()]?)(]?((((](QWWWBBObVW(        
        ´+{[[](?(?)???!)))()(()][{}uhuzzz}[[(}}i!(((??)!(?l???(??i¹:+;+;<l))?!()?)??!Y#$$BBYt),         
       .¹,°;++!!?)!l?i?l!l)?(!([}ubVVYYmaau{}zau}](??!!??!??!)??+,,..-~~~:;>l)?)?l??)hR#@@@X)°          
       )<:;°,-!!?ll?l!llll!???]}{tVMWBWXMOVttbthzz]())()?!???l!!;,´´´´´.,°°~:;+il!?i?uOXB##[i           
       OZOth}!<>li<i>iill?ll!(]}bYZRB$$###BWRZQtuuz{[}[]]])?)?!()!l~.    ´.,¹°°~:>il)}mZWBa[´           
       YMMMMZZZVbz]l+<>+?<?l?]]htQMXWBB$$$###$RYmVYtbh}z[(]]????l)l?)l:.´   ´.,¹°°~;<[tZRQ}-            
       ~)aZRMRRMMMMMZYa{!<>>l[uaVZRXWWWB$$$##$XQQXWMZOmbaaz}]]])!lllili[(~´    ´.,¹°:;)YR}+             
           .~i}YXXRXRRRRRMZtb[?(}tMWWWWB$$$##$XMMW###$WXZOtbhu}](!??!i>i<[b)-     ´-?>¹°i}.             
                ´¹+)bMXXRRXXXXXRQtu}{hYMB$$###BXXB#@@@@@#$XMQVauz}[)l?<<ii>iht<.  ,zu{[i+;+;,           
                      ,~lzOXXWXWXWWWBXQthubVW$BWWB#@@@@@@@##$XZVbz[)l!i<i<!>>>}Za(X$MYaz{(lil<°.        
                           .¹>[mX$$$$$$$$BWRVYmYVX#@@@@@@@@##BWMYh{]?!i<>i>>+><>)bh{R$bzbaz[)?!?;,´     
                                ´,;?aMBB$$B$$$BBWZVOOMM$@@@##BRQYu{]i<><<i><>+>+>>>i]bb..;{bz[(?l>+°.   
                                      .°<}VB$$$$$$$$$WXQOOZZRXRYtz])lli<<<>i>><+<il<th~    -lahz])l>).  
                                           ´-;)bR$$$$$$#$$BBRQQObz{)il<ii<>>>>+>>>l]a]       ´:[buzh]~  
                                                 .~izOB#$##$$$$$BWMQYh[?i><<<><>>>?Yt.          ,>uYu,  
                                                      ´¹+]tX######$$$BBXMOa}?i+><<tY°              .´   
                                                            .~luQ$#$$$$$$BBBWMZY{zVl                    
                                                                 .¹+]mX##$$$$$$BBXh´                    
                                                                       ,~!uZ$#$#$B,                     
                                                                            .¹+[mM                      
                                                                                                    
