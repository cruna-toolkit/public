After successful [installation](install.md) - including the setup of a [cruna.conf](../../cruna.conf.example) - you can run your $CASE as follows:

## 1.) Run CRUNA

`make check`  
`source runtime_env`  
`cp $CAS/parameter.dat $CAS/boundary_conditions.dat bin/`  
`cd ./bin`  
`mpirun -n 16 ./cruna`

## 2.) Postprocess data

### Naming conventions
> Flow field `data_direct_snapshot__01_00001_01_000000.h5`

| name     | block   | proc.   | set      | time step | file type |
| -------- | ------- | ------- | -------- | -------   | -------   |
| data_direct_snapshot | __## | _##### | _## | _###### | .h5 |




