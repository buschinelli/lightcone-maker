def elips_fourier_dist2D(box, spacing):
    ##### Box specifications ######
    Nx, Ny = box.shape[0], box.shape[1]
    dx, dy = spacing[0], spacing[1]
    ##### Values for vec(k) and its components acording to the size and spacing of each real entry

    kx = fftshift(fftfreq(Nx,dx))
    ky = fftshift(fftfreq(Ny,dy))

    ##### Final result matrix 

    ans = np.zeros((Nx,Ny))
    

    ##### Code 

    for i in range(Nx):
        for j in range(Ny):
            ans[i,j] = np.sqrt((kx[i]**2+ky[j]**2))
                


    return ans

def cylindrical_bin(data, spacing, hbins, vbins):
    Nx, Ny, Nz = data.shape[0], data.shape[1], data.shape[2]
    dx, dy, dz = spacing[0], spacing[1], spacing[2]
    spacing2D = np.array([dx, dy])
    
    ans = np.zeros((vbins, hbins))
    kperp_dist = elips_fourier_dist2D(data[:,:,0], spacing2D)
    kpar_max = np.amax(fftshift(fftfreq(Nz,dz)))
    kperp_max = np.amax(kperp_dist)
    delta_kperp = kperp_max/hbins
    delta_kpar = kpar_max/vbins
    kperp_plot = (delta_kperp/2)*np.arange(1,(2*hbins+1),2)
    kpar_plot = (delta_kpar/2)*np.arange(1,(2*vbins+1),2)
    count_bin = np.zeros((vbins,hbins))
    
    for k in range(Nz):
        for i in range(Nx):
            for j in range(Ny):
                kperp = np.sqrt(((i-Nx//2)*(1/(Nx*dx)))**2 + ((j-Ny//2)*(1/(Ny*dy)))**2)
                kpar = (k-Nz//2)*(1/(Nz*dz))
                loc_hbin = int(kperp // delta_kperp)
                loc_vbin = int(kpar // delta_kpar)
                if kperp >= kperp_max:
                    loc_hbin = hbins - 1
                if loc_vbin >= vbins or loc_vbin < -vbins:
                    loc_vbin = vbins - 1
                ans[loc_vbin,loc_hbin] += data[i,j,k]
                count_bin[loc_vbin,loc_hbin] += 1
                
    return ans/count_bin, kpar_plot, kperp_plot

