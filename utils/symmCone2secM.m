function [secM, symmConePara] = symmCone2secM(Emitter)
%symmCone2secM maps the Emitter object (characterized by
%\theta,\phi and rotMobil to orientational second-order moments
mux = sin(Emitter.theta) .* cos(Emitter.phi);
muy = sin(Emitter.theta) .* sin(Emitter.phi);
muz = cos(Emitter.theta);
rotMobility = Emitter.rotMobility;

muxx = rotMobility .* mux.^2 + (1 - rotMobility) / 3;
muyy = rotMobility .* muy.^2 + (1 - rotMobility) / 3;
muzz = rotMobility .* muz.^2 + (1 - rotMobility) / 3;
muxy = rotMobility .* mux .* muy;
muxz = rotMobility .* mux .* muz;
muyz = rotMobility .* muy .* muz;

symmConePara = [mux', muy', rotMobility'];
secM = [muxx', muyy', muzz', muxy', muxz', muyz'];
end