# Institut Fresnel 2001 coherent-data subset

`dielTM_dec8f.exp` is an unchanged local copy of the official Institut Fresnel / IOP 2001 supplementary measurement file that had already been retrieved in this workspace on 2026-09-06. It was copied here on 2026-09-11 because the current IOP page presents a CAPTCHA; no CAPTCHA, registration form, access-control bypass, or external contact was used.

- Object: off-centred homogeneous dielectric circular cylinder (nominal two-dimensional experiment).
- Polarization: TM, electric field parallel to the cylinder axis.
- Frequencies: 1–8 GHz.
- Acquisition: 36 source views; 49 receiver samples per view; 14,112 data rows.
- Format: ASCII comments followed by seven numeric columns: `view receiver frequency_GHz Re(total) Im(total) Re(incident) Im(incident)`. The third column is the physical frequency in GHz; its values happen to be the integers 1 through 8 in this file. The coherent scattered field is `total - incident` under the source paper's `exp(+i omega t)` convention.
- Scientific-use statement: the Institut Fresnel official database page says the first-opus experimental data are free for scientific use and downloadable. The IOP supplementary-data page states that rights remain with the authors unless otherwise specified; this directory does not relicense the data.
- Parser feasibility: high. Whitespace-delimited numeric rows, fixed seven-column layout, finite values, unique `(view, receiver, frequency)` keys. Existing parser examples are available in `research/delegated/a2_measured/ingestion.py` and the public `MPenaR/FresnelDatabaseLoaders` repository.

Current hash:

```text
476cc9d1cfc98797545ab4adf69302dc5aeb45848a24cf8d7b3d222940cc79eb  dielTM_dec8f.exp
```

Current access limits:

- The current IOP data page is bot-gated and was not traversed.
- The Institut Fresnel 3D download route requires account creation and agreement to keep data personal and obtain organizer approval before divulging associated results. No 3D data was acquired or copied for this task.
- The experiment is coherent measured data, but it is a nominally 2D cylinder benchmark. It cannot validate a three-dimensional Gaussian-surface method by itself.

Primary provenance:

- K. Belkebir and M. Saillard, “Special section: Testing inversion algorithms against experimental data,” *Inverse Problems* 17 (2001) 1565–1571, DOI 10.1088/0266-5611/17/6/301.
- Institut Fresnel 3D Database landing page: https://www.fresnel.fr/3Ddatabase/
- IOP supplementary-data page: https://iopscience.iop.org/article/10.1088/0266-5611/17/6/301/data
