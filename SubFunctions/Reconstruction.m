function [Reconstruction_density,Reconstruction_momentumx,Reconstruction_momentumy,Reconstruction_stressx,Reconstruction_stressy,Reconstruction_shear] = Reconstruction(active_elements,mspoints,p_sp,v_ssp,s_sp,GaussCount,x_gauss,y_gauss,x_sp,le,LOCC,Conserve_mass,Conserve_momentum)
ElemsCount = length(active_elements);
Reconstruction_density = zeros(ElemsCount,4);
Reconstruction_momentumx = zeros(ElemsCount,4);
Reconstruction_momentumy = zeros(ElemsCount,4);
Reconstruction_stressx = zeros(ElemsCount,4);
Reconstruction_stressy = zeros(ElemsCount,4);
Reconstruction_shear = zeros(ElemsCount,4);

for i = 1:ElemsCount
    k = active_elements(i);
    particle_numbers = mspoints{k};
    Density = p_sp(particle_numbers);
    Momentumx = p_sp(particle_numbers).*v_ssp(particle_numbers,1);
    Momentumy = p_sp(particle_numbers).*v_ssp(particle_numbers,2);
    Stressx = s_sp(particle_numbers,1);
    Stressy = s_sp(particle_numbers,2);
    Shear = s_sp(particle_numbers,3);
    
[Reconstruction_density(i,:)] = TLS2D(GaussCount,x_gauss(i,:),y_gauss(i,:),length(particle_numbers),x_sp(particle_numbers,:),le,Density,LOCC(k,:),Conserve_mass(i));
[Reconstruction_momentumx(i,:)] = TLS2D(GaussCount,x_gauss(i,:),y_gauss(i,:),length(particle_numbers),x_sp(particle_numbers,:),le,Momentumx,LOCC(k,:),Conserve_momentum(i,1));
[Reconstruction_momentumy(i,:)] = TLS2D(GaussCount,x_gauss(i,:),y_gauss(i,:),length(particle_numbers),x_sp(particle_numbers,:),le,Momentumy,LOCC(k,:),Conserve_momentum(i,2));
[Reconstruction_stressx(i,:)] = TLS2D(GaussCount,x_gauss(i,:),y_gauss(i,:),length(particle_numbers),x_sp(particle_numbers,:),le,Stressx,LOCC(k,:),-999);
[Reconstruction_stressy(i,:)] = TLS2D(GaussCount,x_gauss(i,:),y_gauss(i,:),length(particle_numbers),x_sp(particle_numbers,:),le,Stressy,LOCC(k,:),-999);
[Reconstruction_shear(i,:)] = TLS2D(GaussCount,x_gauss(i,:),y_gauss(i,:),length(particle_numbers),x_sp(particle_numbers,:),le,Shear,LOCC(k,:),-999);
end