function simpplot(p,t,bcol,icol,expr)

dim=size(p,2);
switch dim
 case 2  
  trimesh(t,p(:,1),p(:,2),0*p(:,1),'facecolor','none','edgecolor','k');
  view(2)
  axis equal
  axis off
  
 case 3  
  if nargin>4 & ~isempty(expr)
    tri_full_surf=surftri(p,t);
    
    x = (p(t(:,1),1)+p(t(:,2),1)+p(t(:,3),1)+p(t(:,4),1))/4;
    y = (p(t(:,1),2)+p(t(:,2),2)+p(t(:,3),2)+p(t(:,4),2))/4;
    z = (p(t(:,1),3)+p(t(:,2),3)+p(t(:,3),3)+p(t(:,4),3))/4;
    a = find(eval(expr));
    tnew = t(a,:);
    
    tri_new=surftri(p,tnew);

    % outer boundary
    tri_b = intersect(tri_full_surf,tri_new,'rows');
    
    % inner faces
    tri_i = setdiff(tri_new,tri_b,'rows');
    
    h1=trimesh(tri_b,p(:,1),p(:,2),p(:,3));
    set(h1,'facecolor',bcol,'edgecolor','k');
    hold on;
    h2=trimesh(tri_i,p(:,1),p(:,2),p(:,3));
    set(h2,'facecolor',icol,'edgecolor','k');
    hold off;
    axis equal;
  else
    tri1=surftri(p,t);
    h=trimesh(tri1,p(:,1),p(:,2),p(:,3));
    hold off
    set(h,'facecolor',bcol,'edgecolor','k');
    axis equal
  end
  
 otherwise
  error('Unimplemented dimension.');
end
