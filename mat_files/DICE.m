function dc = DICE(seg1,seg2)

n_vox1=sum(seg1(:)); 
n_vox2=sum(seg2(:));
commarea=sum(seg1(:) & seg2(:)); 
dc = (2*commarea)/(n_vox1+n_vox2);

end
