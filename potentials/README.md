# Si.sw for LAMMPS

Copy `Si.sw` from your LAMMPS installation:

```text
<LAMMPS>/potentials/Si.sw  ->  potentials/Si.sw
```

Example (adjust path):

```powershell
Copy-Item "C:\Program Files\LAMMPS\potentials\Si.sw" potentials\Si.sw
```

Required for V6 TTM-MD pilot (`pair_style sw`).
