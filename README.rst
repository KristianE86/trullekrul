# trullekrul
MATLAB/Octave scripts for anisotropic mesh adaptation, topology optimization and finite element methods

Topology optimization and anisotropic mesh adaptation for optimal heat conduction in 2D and 3D:

.. code:: MATLAB

   top5000(2e-3,5e-3,0.5,false,[],[],0.1,[],1/20,300,1e-3,pi/5,'fig7c');

   top5000(0.4,7e-2,0.1,true,0.75,0.7,0.5,1.5,1/10,283,1e-3,pi/4,'fig9b'):

See the WCSMO12 conference proceeding for a detailed description: https://drive.google.com/file/d/0B9uPpc2f4SZ5TnBjclJFaERSQWs/view?resourcekey=0-oBnY_Ff10qHW27DA4U164w

The same technology can be used for 2D/3D minimum compliance in linear elasticity:

.. code:: MATLAB

   top5001(7.5e2 ,2e-2,0.5 ,0.3,false,2,1, [],1/20,400   ,1e-2,'fig5b',1.025,5,false);

   top5001(2e3   ,2e-2,0.1 ,0.3,2    ,2,0.5,0.25,1/10,400,1e-3,'fig8a',1.025,5,2);

For more details visit the IMR26 paper: https://drive.google.com/open?id=18DRlj6_-MzhOid0BlXC-C-sJvQAm5KnF . Note that the severe objective oscillations reported for 3D in the original IMR26 paper have been fixed.

Finally the technology can be used in Stokes flow with Darcy numbers down to 1e-9, see https://drive.google.com/file/d/1bxHijGOlWnzOBDGfP1Sx5Xk03C9_iy5s/view?usp=sharing

.. code:: MATLAB

    top5002(-1e3, 0.8  ,1e9, [1 1 0.05 0.05       ], 0.01, 141, 'No2',1);

    top5002(-5e3, 0.5,  1e9, [2 3                 ], 0.01, 1e4, 'No9',5);

    top5002(-1e4, 0.1,  1e7, [                    ], 0.01, 1e3, 'No15',3);

Extremely low volume fractions are also trivial to handle.
