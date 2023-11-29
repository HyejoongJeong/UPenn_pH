for(i = 211; i<=225; i++){
j = d2s(i,0);
file_open = "/Volumes/GoogleDrive/Shared drives/Image analysis/Data/EL4_062520/GFP/Position " + j;
file2_open = file_open + "/";
File.openSequence(file2_open);
run("Particle Tracker 2D/3D", "are=No");
run("Particle Tracker 2D/3D", "radius=10 cutoff=0 per/abs=0.10000 link=4 displacement=5 dynamics=Brownian");
file_save ="/Volumes/GoogleDrive/Shared drives/Image analysis/Data/EL4_062520/GFP Trajectories/Position "+j;
file2_save = file_save + ".csv";
saveAs("Results", file2_save);
close("*");
print(i);
}


