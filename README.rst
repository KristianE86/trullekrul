# trullekrul
MATLAB/Octave scripts for anisotropic mesh adaptation, topology optimization and finite element methods

Topology optimization and anisotropic mesh adaptation for optimal heat conduction in 2D and 3D:

.. code:: MATLAB

   top5000(2e-3,5e-3,0.5,false,[],[],0.1,[],1/20,300,1e-3,pi/5,'fig7c')

   top5000(0.4,7e-2,0.1,true,0.75,0.7,0.5,1.5,1/10,283,1e-3,pi/4,'fig9b')

See the WCSMO12 conference proceeding for a detailed description: https://drive.google.com/file/d/0B9uPpc2f4SZ5TnBjclJFaERSQWs

The same technology can be used for 2D/3D minimum compliance in linear elasticity:

.. code:: MATLAB

   top5001(2e-2,1e-2,0.5,0.3,false,2,0.5,[]  ,1/20,200,1e-3,'fig1',1.025,5);

   top5001(4e-2,5e-3,0.1,0.3,1   , 2,0.5,0.25,1/10,400,1e-3,'fig2',1.025,5);

For more details visit the IMR26 paper: http://imr.sandia.gov/papers/imr26/1014_imr26_Jensen.pdf (note that 3D objective oscillations have been fixed)

