function [V,varargout]=variation_operator(mesh,type)
%constructs the variation operator needed for the VB-SCCD algorithm
%
%V=variation_operator(mesh,type)
%
%constructs the variation operator employed in the VB-SCCD source imaging
%algorithm based on the cortical mesh for the case where each triangle 
%corresponds to one grid dipole (type 'face') or the case where each vertex
%corresponds to a grid dipole (type 'vertex').
%
%INPUT: mesh - Matlab structure characterizing the mesh that describes the 
%              triangularized cortical surface
%       optional:
%       type - selects the type of dipole grid:
%               face - each triangle of the mesh corresponds to one dipole
%                      (default)
%               vertex - each vertex of the mesh corresponds to one dipole
%
%OUTPUT: V - E x D matrix of the variation operator (E: number of edges, D:
%            number of grid dipoles)
%
%Hanna Becker, August 2013

%set default type
if nargin==1
    type='face';
end

%Compute number of edges of the mesh
E=mesh.nbfaces*1.5;

%initialization of the variation operator
switch type
    case 'vertex'
        V=sparse(E,length(mesh.v));
    case 'face'
        V=sparse(E,length(mesh.f));
end
adj=zeros(E,2);

%sort indices of vertices for each face
mesh.f=sort(mesh.f,2);

%initialization of edge index
idx=0;
switch type
    case 'face'
        ID=cell(length(mesh.v),1);
        
        %loop of number of faces
        for k=1:mesh.nbfaces
            %for the three pairs of vertices defining an edge of the
            %face...
            
            %first edge
            %determine indices of vertices
            id1=mesh.f(k,1);
            id2=mesh.f(k,2);
            
            %Verify if the edge has been encountered before
            if isempty(ID{id1}) || isempty(find(ID{id1}(1,:)==id2))
                %for a new edge, increase the index for the variation map
                idx=idx+1;
                %store information of edge associated to the vertex with
                %the smallest index
                ID{id1}=[ID{id1},[id2;idx;k;0]];
                %set first element of the variation operator for the edge
                %to 1
                V(idx,k)=1;   
                adj(idx,1)=k;
            else
                %for an already encountered edge, recover the corresponding
                %index of the variation operator from ID and set the second
                %element of the variation operator to -1
                idm=find(ID{id1}(1,:)==id2);
                V(ID{id1}(2,idm),k)=-1;
                adj(ID{id1}(2,idm),2)=k;
                %store information about second edge
                ID{id1}(4,idm)=k;
            end
            
            %second edge
            id1=mesh.f(k,1);
            id2=mesh.f(k,3);
            if isempty(ID{id1}) || isempty(find(ID{id1}(1,:)==id2))
                idx=idx+1;
                ID{id1}=[ID{id1},[id2;idx;k;0]];
                V(idx,k)=1;
                adj(idx,1)=k;
            else
                idm=find(ID{id1}(1,:)==id2);
                V(ID{id1}(2,idm),k)=-1;
                adj(ID{id1}(2,idm),2)=k;
                ID{id1}(4,idm)=k;
            end
            
            %third edge
            id1=mesh.f(k,2);
            id2=mesh.f(k,3);
            if isempty(ID{id1}) || isempty(find(ID{id1}(1,:)==id2))
                idx=idx+1;
                ID{id1}=[ID{id1},[id2;idx;k;0]];
                V(idx,k)=1;    
                adj(idx,1)=k;
            else
                idm=find(ID{id1}(1,:)==id2);
                V(ID{id1}(2,idm),k)=-1;
                adj(ID{id1}(2,idm),2)=k;
                ID{id1}(4,idm)=k;
            end
        end
        varargout{1}=adj;               
        
    case 'vertex'        
        ID=cell(length(mesh.v),1);
        for k=1:mesh.nbfaces
            
            id1=mesh.f(k,1);
            id2=mesh.f(k,2);
            if isempty(find(ID{id1}==id2))
                idx=idx+1;
                ID{id1}=[ID{id1},id2];
                V(idx,id1)=1;
                V(idx,id2)=-1;
            end
            
            id1=mesh.f(k,1);
            id2=mesh.f(k,3);
            if isempty(find(ID{id1}==id2))
                idx=idx+1;
                ID{id1}=[ID{id1},id2];
                V(idx,id1)=1;
                V(idx,id2)=-1;
            end
            
            id1=mesh.f(k,2);
            id2=mesh.f(k,3);
            if isempty(find(ID{id1}==id2))
                idx=idx+1;
                ID{id1}=[ID{id1},id2];
                V(idx,id1)=1;
                V(idx,id2)=-1;
            end
        end        
end