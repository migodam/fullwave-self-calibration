# Extension note: not part of the Q5 manuscript or its evidence

Heartbeat/respiration separation can be posed as a different task-conditioned inverse problem. If a physical dynamic scattering model is available, write y(t)=g(t)F(q_resp(t),q_heart(t),s(t),epsilon)+n(t). Material and electronic calibration are not automatically substitutes for identifying the two motions. Respiration harmonics, common body displacement and multipath may overlap the nominal cardiac band.

The transferable idea is to define the two desired motion tasks and the admissible full trajectories, then ask whether measurement-compatible trajectory pairs can differ by more than the task tolerance. A frequency filter is sufficient only under proved nonoverlap/aliasing assumptions. A physical motion derivative may be used for a local diagnostic; an arbitrary nuisance vector cannot be labelled respiration. A shared electronic reference helps only the shared electronic nuisance. Additional look directions, polarizations, temporal priors or reference motion measurements need their own acquisition/noise/forward-error model.

No physiological data, human subject, clinical validation, new heartbeat algorithm or medical performance claim is part of V2. There is no demonstrated engineering transfer yet. The code and proofs remain restricted to the declared electromagnetic sphere problems.
