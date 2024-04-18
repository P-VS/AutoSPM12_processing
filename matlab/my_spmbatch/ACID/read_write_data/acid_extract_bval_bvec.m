function [bval,bvec] = acid_extract_bval_bvec(P, filename, p_out)
% B.Fricke 17/12/2021

P = unique(P,'rows');

try PJ = acid_get_metadata(P);

    [bval,bvec] = read_bval_bvec_json(PJ,filename,p_out);
    normgrad = sqrt(sum(bvec .* bvec, 1));
    tolgrad    = 0.01; % tolerance for normalizing gradient directions
    MSKgrad = find(normgrad > 0);

    if(~isempty(find(normgrad(MSKgrad) ~= 1,1)))
        if(~isempty(find(abs(normgrad(MSKgrad)-1) > tolgrad,1))) % SM: if norm of gradients is within the tolerance message won't be visible / BF: And no changes on the b-vectors
            warning('Diffusion gradients are now normalised!');
            bvec(:, MSKgrad) = bsxfun(@rdivide, bvec(:, MSKgrad), normgrad(MSKgrad));
        end
    end
    clear normgrad;

    if size(P,1) == size(bval,2) && size(bval,2) == size(bvec,2)
        V4D = spm_file_merge(P,[p_out filesep filename '.nii']);

        

        save([p_out filesep filename '_bval-bvec.mat'],'bval','bvec')
        Mtmp1 = [p_out filesep filename '.bval'];
        Mtmp2 = [p_out filesep filename '.bvec'];
        dlmwrite(Mtmp1, bval);
        dlmwrite(Mtmp2, bvec);
    elseif size(P,1) ~= size(bval,2)
        error('Number of b-values and images are not the same!')
    elseif size(P,1) ~= size(bvec,2)
        error('Number of b-vectors and images are not the same!')
    else
        error('Number of b-vectors and b-values are not the same!')
    end

catch
    error('No bval/bvec file is created!')
end
end

function [bval,bvec] = read_bval_bvec_json(PJ,fname,p_out)
% Function to extract b-value and b-vectors from json header.
% S.Mohammadi 05/11/2021
% B.Fricke 17/12/2021

bvec = zeros(3,size(PJ,2));
bval = zeros(1,size(PJ,2));

for inx = 1:size(PJ,2)

    if isfield(PJ{inx}.acqpar,'ImplementationVersionName') && strcmp(PJ{inx}.acqpar.ImplementationVersionName,'9_R292a_M4_Patch')

        try bval(inx) = PJ{inx}.acqpar.Private_0043_1039(1,1);
            if ~strcmp(PJ{inx}.acqpar.Private_0019_10bb,'NONE')
                bvec(1,inx) = PJ{inx}.acqpar.Private_0019_10bb;
                bvec(2,inx) = PJ{inx}.acqpar.Private_0019_10bc;
                bvec(3,inx) = PJ{inx}.acqpar.Private_0019_10bd;
            else
                bvec(:,inx) = zeros(3,1);
            end
        catch
        end

    elseif isfield(PJ{inx}.acqpar,'SoftwareVersions') && (strcmp(PJ{inx}.acqpar.SoftwareVersions,'syngo MR B17' ) ||  strcmp(PJ{inx}.acqpar.SoftwareVersions,'syngo MR E12'))

        try bval(inx) = PJ{inx}.acqpar.CSAImageHeaderInfo.B_value(1,1);
            if isfield(PJ{inx}.acqpar.CSAImageHeaderInfo,'DiffusionGradientDirection')
                bvec(1,inx) = PJ{inx}.acqpar.CSAImageHeaderInfo.DiffusionGradientDirection(1,1);
                bvec(2,inx) = PJ{inx}.acqpar.CSAImageHeaderInfo.DiffusionGradientDirection(2,1);
                bvec(3,inx) = PJ{inx}.acqpar.CSAImageHeaderInfo.DiffusionGradientDirection(3,1);
            else
                bvec(:,inx) = zeros(3,1);
            end
        catch
        end
    elseif isfield(PJ{inx}.acqpar,'ManufacturerModelName') && strcmp(PJ{inx}.acqpar.ManufacturerModelName,'MAGNETOM Prisma ')

        try bval(inx) = PJ{inx}.acqpar.PerFrameFunctionalGroupsSequence(1).MRDiffusionSequence.DiffusionBValue(1,1);
            if isfield(PJ{inx}.acqpar.PerFrameFunctionalGroupsSequence(1).MRDiffusionSequence,'DiffusionGradientDirectionSequence')
                bvec(1,inx) = PJ{inx}.acqpar.PerFrameFunctionalGroupsSequence(1).MRDiffusionSequence.DiffusionGradientDirectionSequence.DiffusionGradientOrientation(1,1);
                bvec(2,inx) = PJ{inx}.acqpar.PerFrameFunctionalGroupsSequence(1).MRDiffusionSequence.DiffusionGradientDirectionSequence.DiffusionGradientOrientation(2,1);
                bvec(3,inx) = PJ{inx}.acqpar.PerFrameFunctionalGroupsSequence(1).MRDiffusionSequence.DiffusionGradientDirectionSequence.DiffusionGradientOrientation(3,1);
            else
                bvec(:,inx) = zeros(3,1);
            end
        catch
        end
    elseif isfield(PJ{inx}.acqpar,'ImplementationVersionName') && (strcmp(PJ{inx}.acqpar.ImplementationVersionName,'SYNGO_MR_XA61A') ||  strcmp(PJ{inx}.acqpar.ImplementationVersionName,'SYNGO_MR_XA60A') )

        try bval(inx) = PJ{inx}.acqpar.PerFrameFunctionalGroupsSequence(1).MRDiffusionSequence.DiffusionBValue(1,1);
            if isfield(PJ{inx}.acqpar.PerFrameFunctionalGroupsSequence(1).MRDiffusionSequence,'DiffusionGradientDirectionSequence')
                bvec(1,inx) = PJ{inx}.acqpar.PerFrameFunctionalGroupsSequence(1).MRDiffusionSequence.DiffusionGradientDirectionSequence.DiffusionGradientOrientation(1,1);
                bvec(2,inx) = PJ{inx}.acqpar.PerFrameFunctionalGroupsSequence(1).MRDiffusionSequence.DiffusionGradientDirectionSequence.DiffusionGradientOrientation(2,1);
                bvec(3,inx) = PJ{inx}.acqpar.PerFrameFunctionalGroupsSequence(1).MRDiffusionSequence.DiffusionGradientDirectionSequence.DiffusionGradientOrientation(3,1);
            else
                bvec(:,inx) = zeros(3,1);
            end
        catch
        end
    else
        try bval(inx) = PJ{inx}.acqpar.Private_0019_100c;

            if ~strcmp(PJ{inx}.acqpar.Private_0019_100d,'NONE')
                bvec(:,inx) = PJ{inx}.acqpar.Private_0019_100e;
            else
                bvec(:,inx) = zeros(3,1);
            end

        catch

            try bval(inx) = PJ{inx}.acqpar.DiffusionBValue;
                if ~strcmp(PJ{inx}.acqpar.Private_0019_10bb,'NONE')
                    bvec(1,inx) = PJ{inx}.acqpar.Private_0019_10bb;
                    bvec(2,inx) = PJ{inx}.acqpar.Private_0019_10bc;
                    bvec(3,inx) = PJ{inx}.acqpar.Private_0019_10bd;
                else
                    bvec(:,inx) = zeros(3,1);
                end
            catch
            end
        end
    end
end



% Rotating bvecs in dependency of patient orientation
try rotation_matrix = [PJ{1}.acqpar.ImageOrientationPatient(1:3)';PJ{1}.acqpar.ImageOrientationPatient(4:6)';0,0,1];

    rotation_matrix(3,:)= cross(rotation_matrix(1,:)',rotation_matrix(2,:)'); % Or json: "SliceNormalVector"

    bvec = rotation_matrix*bvec;
catch

    warning('Rotation parameters for bvecs not found! bvecs might not be in patient orientation system!')

end

refNtmp = round(bval/100)*100;
[aa,~,~] = unique(refNtmp);
if aa(1,1) == 0
    bb = histc(refNtmp,aa(1,2:end));
else
    bb = histc(refNtmp,aa);
end

bb_count = histc(refNtmp,aa);

if min(bb) < 3
    warning('One diffusion shell contains less then three diffusion directions!')
end

parameterstruct = PJ{1}.acqpar;
textfile=fopen([p_out filesep fname '-Parameter.txt'],'w');
fprintf(textfile,'%s\n','General parameters:');

try fprintf(textfile,'%s\n','  ');
catch
end
try fprintf(textfile,'%s\n',['Manufacturer: ' parameterstruct.Manufacturer]);
catch
end
try fprintf(textfile,'%s\n',['Model name: ' parameterstruct.ManufacturerModelName]);
catch
end
try fprintf(textfile,'%s\n',['Institution name: ' parameterstruct.InstitutionName]);
catch
end
try fprintf(textfile,'%s\n',['Series description: ' parameterstruct.SeriesDescription]);
catch
end

fprintf(textfile,'%s\n','  ');
fprintf(textfile,'%s\n','  ');
fprintf(textfile,'%s\n','  ');

fprintf(textfile,'%s\n','Patient informations:');
fprintf(textfile,'%s\n','  ');

try fprintf(textfile,'%s\n',['Patient ID: ' parameterstruct.PatientID]);
catch
end
try fprintf(textfile,'%s\n',['Patient sex: ' parameterstruct.PatientSex]);
catch
end
try fprintf(textfile,'%s\n',['Patient age: ' num2str(parameterstruct.PatientAge)]);
catch
end
try fprintf(textfile,'%s\n',['Patient size: ' num2str(parameterstruct.PatientSize)]);
catch
end
try fprintf(textfile,'%s\n',['Patient weight: ' num2str(parameterstruct.PatientWeight)]);
catch
end

fprintf(textfile,'%s\n','  ');
fprintf(textfile,'%s\n','  ');
fprintf(textfile,'%s\n','  ');

fprintf(textfile,'%s\n','Acquisition parameters:');
fprintf(textfile,'%s\n','  ');

try fprintf(textfile,'%s\n',['Acquisition time: ' num2str(parameterstruct.AcquisitionTime)]);
catch
end
try fprintf(textfile,'%s\n',['MR acquisition type: ' parameterstruct.MRAcquisitionType]);
catch
end
try fprintf(textfile,'%s\n',['Sequence name: ' parameterstruct.SequenceName]);
catch
end
try fprintf(textfile,'%s\n',['Field of View: ' num2str(parameterstruct.Private_0051_100c)]);
catch
end
try fprintf(textfile,'%s\n',['Slice thickness: ' num2str(parameterstruct.SliceThickness)]);
catch
end
try fprintf(textfile,'%s\n',['Repetition Time (TR): ' num2str(parameterstruct.RepetitionTime)]);
catch
end
try fprintf(textfile,'%s\n',['Repetition Time (TR): ' num2str(parameterstruct.SharedFunctionalGroupsSequence.MRTimingAndRelatedParametersSequence.RepetitionTime)]);
catch
end
try fprintf(textfile,'%s\n',['Echo Time (TE): ' num2str(parameterstruct.EchoTime)]);
catch
end
try fprintf(textfile,'%s\n',['Effective Echo Time (TE): ' num2str(parameterstruct.PerFrameFunctionalGroupsSequence(1).MREchoSequence.EffectiveEchoTime)]);
catch
end
try fprintf(textfile,'%s\n',['Effective Train length: ' num2str(parameterstruct.SharedFunctionalGroupsSequence.MRTimingAndRelatedParametersSequence.EchoTrainLength)]);
catch
end
try fprintf(textfile,'%s\n',['Number of averages: ' num2str(parameterstruct.NumberOfAverages)]);
catch
end
try fprintf(textfile,'%s\n',['Echo numbers: ' num2str(parameterstruct.EchoNumbers)]);
catch
end
try fprintf(textfile,'%s\n',['Magnetic field strength: ' num2str(parameterstruct.MagneticFieldStrength)]);
catch
end
try fprintf(textfile,'%s\n',['Spacing between slices: ' num2str(parameterstruct.SpacingBetweenSlices)]);
catch
end
try fprintf(textfile,'%s\n',['Number of phase encoding steps: ' num2str(parameterstruct.NumberOfPhaseEncodingSteps)]);
catch
end
try fprintf(textfile,'%s\n',['Echo train length: ' num2str(parameterstruct.EchoTrainLength)]);
catch
end
try fprintf(textfile,'%s\n',['Percent sampling: ' num2str(parameterstruct.PercentSampling)]);
catch
end
try fprintf(textfile,'%s\n',['Percent phase field of view: ' num2str(parameterstruct.PercentPhaseFieldOfView)]);
catch
end
try fprintf(textfile,'%s\n',['Protocol name: ' parameterstruct.ProtocolName]);
catch
end
try fprintf(textfile,'%s\n',['Transmit coil name: ' parameterstruct.TransmitCoilName]);
catch
end
try fprintf(textfile,'%s\n',['Acquisition matrix: ' num2str(parameterstruct.AcquisitionMatrix(:,1)')]);
catch
end
try fprintf(textfile,'%s\n',['In plane phase encoding direction: ' parameterstruct.InPlanePhaseEncodingDirection]);
catch
end
try fprintf(textfile,'%s\n',['Flip angle: ' num2str(parameterstruct.FlipAngle)]);
catch
end
try fprintf(textfile,'%s\n',['Flip angle: ' num2str(parameterstruct.SharedFunctionalGroupsSequence.MRTimingAndRelatedParametersSequence.FlipAngle)]);
catch
end
try fprintf(textfile,'%s\n',['SAR: ' num2str(parameterstruct.SAR)]);
catch
end

fprintf(textfile,'%s\n','  ');
fprintf(textfile,'%s\n','  ');
fprintf(textfile,'%s\n','  ');

try fprintf(textfile,'%s\n',['Measured shells:                    ' num2str(aa)]);
catch
end

try fprintf(textfile,'%s\n',['Number of b-vectors for each shell: ' num2str(bb_count)]);
catch
end

fclose(textfile);
end