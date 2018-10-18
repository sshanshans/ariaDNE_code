function Write(G,filename,format,options)
% Color added to .off files by Julie Winchester (julie.winchester@duke.edu)

options.pointCloud = getoptions(options, 'pointCloud', 0);
options.color = getoptions(options, 'color', []);

switch format
    case 'off'
        fid = fopen(filename,'wt');
        if( fid==-1 )
            error('Can''t open the file.');
        end
        
        % header
        fprintf(fid, 'OFF\n');
        if options.pointCloud==0
            fprintf(fid, '%d %d 0\n', length(G.V), length(G.F));
        else
            fprintf(fid, '%d 0 0\n', length(G.V));
        end
        
        % write the points & faces
        fprintf(fid, '%f %f %f\n', G.V);
        if options.pointCloud==0
            if ~isempty(options.color)
                fc = vertcat(G.F - 1, options.color);
                fprintf(fid, ['3 %d %d %d %f %f %f %f\n'], fc);
            else
                fprintf(fid, ['3 %d %d %d\n'], G.F - 1);
            end
        end
        
        fclose(fid);
    case 'ply'
        fid = fopen(filename,'wt');
        if( fid==-1 )
            error('Can''t open the file.');
        end

        fprintf(fid, 'ply\nformat ascii 1.0\n');
        fprintf(fid, 'element vertex %d\n', length(G.V));
        fprintf(fid, 'property float x\nproperty float y\nproperty float z\n');
        if ~isempty(options.color)
            fprintf(fid, 'property uchar red\nproperty uchar green\nproperty uchar blue\nproperty uchar alpha\n');
        end
        if options.pointCloud==0
            fprintf(fid, 'element face %d\n', length(G.F));
            fprintf(fid, 'property list uchar int vertex_indices\n');
        end
        fprintf(fid, 'end_header\n');

        % write the points & faces
        if ~isempty(options.color)
            c = repmat(options.color(:,1), [1, length(G.V) ]) * 255;
            vc = vertcat(G.V, c);
            fprintf(fid, '%f %f %f %d %d %d %d\n', vc);
        else
            fprintf(fid, '%f %f %f\n', G.V);
        end
            
        if options.pointCloud==0
            fprintf(fid, ['3 %d %d %d\n'], G.F - 1);
        end

    case 'obj'
        fid = fopen(filename,'wt');
        if( fid==-1 )
            error('Can''t open the file.');
        end
        
        % vertex coordinates
        fprintf(fid, 'v %f %f %f\n', G.V);
        
        % Texture coordinates
        if isfield(options, 'Texture')
            fprintf(fid, 'vt %f %f\n', options.Texture.Coordinates(1:2,:));
            fprintf(fid, 'f %d/%d %d/%d %d/%d\n', kron(G.F',[1,1])');
        else
            fprintf(fid, 'f %d %d %d\n', G.F');
        end
        fclose(fid);
end