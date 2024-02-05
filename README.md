# master-project-final

## Prerequisites

Please check how to build ns-3 Projects and Prerequisites here:

> https://github.com/nsnam/ns-3-dev-git/blob/master/doc/installation/source/quick-start.rst

## Build

**1) Clone Repository**
-----------------------

```console
git clone https://github.com/marbuntu/master-project-final.git
cd master-project-final
```

**2) Configure**
------------------------
```
./waf configure
```

**3) Build**
-------------
```
./waf build
```

**4) Run**
-----------
ns-3 NTN Simulation:
```
./waf --run mpj-simulation
```

Constellation Helper Test
```
./waf --run usr-constellation-test
```



## References
---------------
Original ns-3 Repo

> https://github.com/nsnam/ns-3-dev-git


sns-3 Satellite

> https://github.com/sns3/sns3-satellite

ns-3 satellite

> https://gitlab.inesctec.pt/pmms/ns3-satellite

n3-3 NTN

> https://gitlab.com/mattiasandri/ns-3-ntn