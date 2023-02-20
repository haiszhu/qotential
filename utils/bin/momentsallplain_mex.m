function m_all=momentsallplain_mex(r0,r,order,np,p)


m = numel(r0(1,:));
n = numel(r(1,:));
dim1 = n*(2*order+2);

mex_id_ = 'momentsallplain(i int[x], i double[xx], i int[x], i double[xx], i int[x], i int[x], i int[x], i int[x], o double[xx])';
[m_all] = specialquad(mex_id_, m, r0, n, r, order, np, p, dim1, 1, 3, m, 1, 3, n, 1, 1, 1, 1, dim1, m);

m_all = reshape(m_all,n,2*order+2,m);

end

% @function [n_all,m_all,l_all]=momentsall_mex(r0,r,order,sqn_flag,r_root,t_root,tj,wj,D,w_bclag,np,p)
% 
% nout = double(nargout);
% 
% m = numel(r0(1,:));
% n = numel(r(1,:));
% sqn_flag = double(sqn_flag);
% np3 = np*3;
% r_root = reshape(r_root,np3,m);
% dim1 = n*(2*order+2);
% 
% if nout < 3
% # FORTRAN momentsall(int[1] m, double[3,m] r0, int[1] n, double[3,n] r, int[1] order, int[np,m] sqn_flag, double[np3,m] r_root, double[np,m] t_root, double[p] tj, double[p] wj, double[p,p] D, double[p] w_bclag, int[1] np, int[1] p, int[1] dim1, output double[dim1,m] n_all, output double[dim1,m] m_all);
% 
%   m_all = reshape(m_all,n,2*order+2,m);
%   n_all = reshape(n_all,n,2*order+2,m);
% else
% # FORTRAN moments2all(int[1] m, double[3,m] r0, int[1] n, double[3,n] r, int[1] order, int[np,m] sqn_flag, double[np3,m] r_root, double[np,m] t_root, double[p] tj, double[p] wj, double[p,p] D, double[p] w_bclag, int[1] np, int[1] p, int[1] dim1, output double[dim1,m] n_all, output double[dim1,m] m_all, output double[dim1,m] l_all);
% 
%   l_all = reshape(l_all,n,2*order+2,m);
%   m_all = reshape(m_all,n,2*order+2,m);
%   n_all = reshape(n_all,n,2*order+2,m);
% end
% 
% % keyboard
% 
% end

% @function mk_j=momentsup_mex(tgl,wgl,r,r0,t_root,r_root,w_bclag,D,order)
% 
% p = numel(tgl);
% morderp1 = 2*order + 2;
% # FORTRAN momentsup(int[1] p, double[p] tgl, double[p] wgl, double[3,p] r, double[3] r0, double[1] t_root, double[3] r_root, double[p] w_bclag, double[p,p] D, int[1] order, output double[p,morderp1] nk_j, output double[p,morderp1] mk_j);
% 
% if 0 % separate routine...
% % compute length of pan_t_end
% # FORTRAN lenofpantend(int[1] p, double[p] tgl, double[p] wgl, double[3,p] r, double[3] r0, double[1] t_root, double[3] r_root, double[p] w_bclag, double[p,p] D, output double[1] sqn_dist, output double[1] pan_len, output int[1] len, output int[1] lenl, output int[1] lenr, output double[3,p] rp, output double[p] w);
% n = p*(len-1);
% # FORTRAN momentsup0(int[1] p, double[p] tgl, double[p] wgl, double[3,p] r, double[3] r0, double[1] t_root, double[3] r_root, double[p] w_bclag, double[p,p] D, int[1] order, double[1] sqn_dist, double[1] pan_len, int[1] len, int[1] lenl, int[1] lenr, double[3,p] rp, double[p] w, int[1] n, output double[p,morderp1] nk_j, output double[p,morderp1] mk_j, output double[1] test);
% end
% 
% end


% @function [r_up,Br]=line3adaptivernearr0_mex(tgl,wgl,r,r0,t_root,r_root,w_bclag,D)
% 
% p = numel(tgl);
% 
% % compute length of pan_t_end
% # FORTRAN lenofpantend(int[1] p, double[p] tgl, double[p] wgl, double[3,p] r, double[3] r0, double[1] t_root, double[3] r_root, double[p] w_bclag, double[p,p] D, output double[1] sqn_dist, output double[1] pan_len, output int[1] len, output int[1] lenl, output int[1] lenr, output double[3,p] rp, output double[p] w);
% 
% % pan_t_end
% # FORTRAN addextrapanel(double[1] sqn_dist, double[1] pan_len, double[1] t_root, int[1] len, int[1] lenl, int[1] lenr,  output double[len] pan_t_end, output int[1] npl);
% 
% n = p*(len-1);
% # FORTRAN line3adaptivernearr0(int[1] p, double[p] tgl, double[p] wgl, double[3,p] r, double[p] w_bclag, int[1] len, double[3,p] rp, double[p] w, double[len] pan_t_end, int[1] n, output double[3,n] r_up, output double[3,n] rp_up, output double[p,n] Br, output double[n] test);
% 
% end

% @function [sqn_dist,pan_len,len,lenl,lenr,rp,w]=lenofpantend_mex(tgl,wgl,r,r0,t_root,r_root,w_bclag,D)
% 
% m = numel(tgl);
% 
% # FORTRAN lenofpantend(int[1] m, double[m] tgl, double[m] wgl, double[3,m] r, double[3] r0, double[1] t_root, double[3] r_root, double[m] w_bclag, double[m,m] D, output double[1] sqn_dist, output double[1] pan_len, output int[1] len, output int[1] lenl, output int[1] lenr, output double[3,m] rp, output double[m] w);
% 
% end

% @function [pan_t_end, npl]=addextrapanel_mex(sqn_dist, pan_len, t_root)
% 
% % compute length of pan_t_end
% factor = 1.5;
% if t_root >= 1 % panel on the left side of the target
%   pan_len_l = factor*pan_len;
%   ml = floor(log(pan_len_l/sqn_dist)/log(factor));
%   ml = max(ml,0);
%   mr = 1; % fake, should be 0, no use
%   m = ml+2;
% elseif t_root <= -1 % panel on the right side of the target
%   pan_len_r = factor*pan_len;
%   mr = floor(log(pan_len_r/sqn_dist)/log(factor));
%   mr = max(mr,0);
%   ml = 1; % fake, should be 0, no use
%   m = mr+2;
% else  % subdivide both sides
%   pan_len_l = factor*(t_root-(-1))/2*pan_len; % roughtly in physical space
%   ml = floor(log(pan_len_l/sqn_dist)/log(factor)); ml = max(ml,0);
%   pan_len_r = factor*(1-t_root)/2*pan_len;
%   mr = floor(log(pan_len_r/sqn_dist)/log(factor)); mr = max(mr,0);
%   m = ml+mr+3;
% end
% 
% # FORTRAN addextrapanel(double[1] sqn_dist, double[1] pan_len, double[1] t_root, int[1] m, int[1] ml, int[1] mr,  output double[m] pan_t_end, output int[1] npl);
% 
% % do any MATLAB postprocessing here
% pan_t_end = pan_t_end';
% end

% @function B=bclaginterpmatrix_mex(x, xi, w)
% n = numel(x);
% m = numel(xi);
% # FORTRAN bclaginterpmatrix(double[n] x, int n, double[m] xi, int m, double[n] w, output double[m,n] B);
% % do any MATLAB postprocessing here
% end

