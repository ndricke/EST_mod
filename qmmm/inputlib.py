#a bunch of variables for writing QM/MM input files in Q-Chem

charge_tran = {'a2':-2,'a1':-1,'a0':0,'c1':1,'c2':2,'c3':3}
#revd=dict([reversed(i) for i in charge_tran.items()])
#charge_tran.update(revd)

qmmm_rem1 = \
"""
$rem
method                    b3lyp
basis                     6-31G*
qm_mm_interface           janus
force_field               oplsaa
user_connect              true
jobtype                   aimd
scf_max_cycles            2000
scf_algorithm             DIIS_GDM
time_step                 20
aimd_steps                2
aimd_init_veloc           thermal
aimd_thermostat           nose_hoover
nose_hoover_length        3
nose_hoover_timescale     10
aimd_temp                 300
"""

##Also removed these; didn't seem to be necessary for now
#molden_format           true
#print_orbitals          true

qmmm_rem2 = \
"""
$rem
method                    b3lyp
basis                     6-31G*
qm_mm_interface           janus
force_field               oplsaa
user_connect              true
scf_guess                 read
scf_max_cycles            2000
scf_algorithm             DIIS_GDM
gaussian_blur             true
gauss_blur_width          15000
skip_scfman               true
jobtype                   aimd
time_step                 20
aimd_steps                100
aimd_init_veloc           thermal
aimd_temp                 300
aimd_thermostat           nose_hoover
nose_hoover_length        3
nose_hoover_timescale     10
"""

qmmm_remvel = \
"""
$rem
method                    b3lyp
basis                     6-31G*
qm_mm_interface           janus
force_field               oplsaa
user_connect              true
scf_max_cycles            2000
scf_algorithm             DIIS_GDM
gaussian_blur             true
gauss_blur_width          15000
jobtype                   aimd
time_step                 1
aimd_steps                1
aimd_init_veloc           thermal
aimd_temp                 300
aimd_thermostat           nose_hoover
nose_hoover_length        3
nose_hoover_timescale     10
"""

forceman = \
"""
$forceman
ewald
alpha 0.2 0.01
box_length 32.0 32.0 32.0
$end
"""

Fe_ff = \
"""
$force_field_params
NumAtomTypes 1
AtomType -1  1.45259  2.93  0.0
$end
"""

read = \
"""
@@@

$molecule
read
$end
"""
