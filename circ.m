classdef circ
    
    properties(Constant)
        colours = {'b','r','g','k','c','m','y'};
    end
    
    properties
        componentList
        voltageNodes
        numComp
        numNodes
        numEqn
        Gmat
        Cmat
        Fmat
        Vmat
        solution
        status
        solvedValuesNames
        inputComponentName
        inputComponentCompNum
        inputComponentEqnNum
        sweepValues
        numSweepValues
        acSweepInput
    end
    
    methods(Static)
        function [out] = makeLength(in,L)
            out = [in,blanks(L-length(in))];
        end
        
        function [tf] = checkInt(val,err)
            if nargin == 1
                err = 0.000001;
            end
            if abs(floor(val) - val) <= err
                tf = true;
            else
                tf = false;
            end
        end
        
        function [out] = formatVal(val,L)
            if nargin == 1
                L = -1;
            end
            if circ.checkInt(val)
                out = sprintf('%d',val);
            else
                if floor(log10(val)) < -2
                    out = sprintf('%.e',val);
                else
                    out = sprintf('%.2f',val);
                end
            end
            out = circ.makeLength(out,L);
        end
        
        function [out] = printVal(val)
            val = abs(val);
            if val == 1
                out = '';
                return;
            end
            out = [circ.formatVal(val),'*'];
        end
        
        function [side] = addToSideIn(side,text)
            if isempty(side)
                side = text;
            else
                side = [side,' + ',text];
            end
        end
        
        function [ls,rs] = addToSide(ls,rs,text,val)
            s = sign(val);
            if s == -1
                ls = circ.addToSideIn(ls,text);
            elseif s == 1
                rs = circ.addToSideIn(rs,text);
            end
        end
        
        function [f,S] = timeToFreq(t,V)
            N = length(V);
            S = abs(fft(V));
            S = S/N;
            S(2:N/2) = S(2:N/2) + flip(S((N/2+2):N));
            S = S(1:N/2+1);
            f = (0:N/2)/(max(t)-min(t));
        end
        
        function [lineSpec] = getLineSpec(n)
            lineSpec = sprintf('%s.-',circ.colours{n});
        end
    end
    
    methods(Access = private)
        
        %% Setup methods
        function obj = performKCL(obj)
            tmp = cell(1,obj.numNodes);
            for nn = 1:obj.numComp
                currComp = obj.componentList(nn);
                if currComp.n1 ~= 0
                    tmp{currComp.n1} = [tmp{currComp.n1},currComp];
                end
                if currComp.n2 ~= 0
                    tmp{currComp.n2} = [tmp{currComp.n2},currComp];
                end
                obj.solvedValuesNames{nn} = sprintf('I(%s)',currComp.name);
            end
            for nn = 1:obj.numNodes
                currNodeComps = tmp{nn};
                obj.solvedValuesNames{nn + obj.numComp} = sprintf('V(%d)',nn);
                for mm = 1:length(currNodeComps)
                    currComp = currNodeComps(mm);
                    s = 1;
                    if currComp.n1 == nn
                        s = -1;
                    end
                    obj.Gmat(nn,currComp.num) = s;
                end
            end
        end
        
        function obj = analyzeComps(obj)
            for currComp = obj.componentList
                currEqn = currComp.num + obj.numNodes;
                currType = currComp.getType();
                if currType == 'v' || currType == 'r' || currType == 'l' || currType == 'b'
                    if currComp.n1 ~= 0
                        obj.Gmat(currEqn,currComp.n1 + obj.numComp) = -1;
                    end
                    if currComp.n2 ~= 0
                        obj.Gmat(currEqn,currComp.n2 + obj.numComp) = 1;
                    end
                elseif currType == 'c' || currType == 'i'
                    obj.Gmat(currEqn,currComp.num) = 1;
                end
                obj = obj.assignVal(currComp);
            end
        end
        
        function [obj] = assignVal(obj,currComp)
            currEqn = currComp.num + obj.numNodes;
            switch currComp.getType()
                case 'v'
                    obj.Fmat(currEqn) = currComp.val;
                case 'r'
                    obj.Gmat(currEqn,currComp.num) = currComp.val;
                case 'l'
                    obj.Cmat(currEqn,currComp.num) = currComp.val;
                case 'c'
                    if currComp.n1 ~= 0
                        obj.Cmat(currEqn,currComp.n1 +obj.numComp) = -currComp.val;
                    end
                    if currComp.n2 ~= 0
                        obj.Cmat(currEqn,currComp.n2 + obj.numComp) = currComp.val;
                    end
                case 'b'
                    obj.Gmat(currEqn,obj.getCompNum(currComp.rel)) = -currComp.val;
                case 'i'
                    obj.Fmat(currEqn) = currComp.val;
            end
        end
        
        %% Solving methods
        
        function [obj] = setSweepValues(obj,sweepValues)
            obj.sweepValues = sweepValues;
            obj.numSweepValues = length(obj.sweepValues);
            obj.solution = zeros(obj.numEqn,obj.numSweepValues);
        end
        
        function [obj] = solve(obj,Vin,w)
            obj = setupIfNotAlready(obj);
            if nargin >= 2
                obj.Fmat(obj.inputComponentEqnNum) = Vin;
            end
            if nargin == 3
                obj.Vmat = mldivide(obj.Gmat + w*1j*obj.Cmat,obj.Fmat);
            elseif nargin == 2
                obj.Vmat = mldivide(obj.Gmat,obj.Fmat);
            else
                obj.Vmat = mldivide(obj.Gmat,obj.Fmat);
            end
        end
        
        %% Output methods
        
        function [] = outputOpPoint(obj,compIndex)
            for n = 1:length(compIndex)
                if imag(obj.solution(n)) == 0
                    fprintf('%s = %f\n',obj.solvedValuesNames{compIndex(n)},obj.solution(compIndex(n)));
                else
                    fprintf('%s = %%f+%fj\n',obj.solvedValuesNames{compIndex(n)},real(obj.solution(compIndex(n))),imag(obj.solution(compIndex(n))));
                end
            end
            fprintf('\n');
        end
        
        function [f] = outputDC(obj,compIndex)
            f = figure;
            for n = 1:length(compIndex)
                lineSpec = circ.getLineSpec(n);
                subplot(1,1,1);
                hold on;
                if obj.isVarVoltage(compIndex(n))
                    yyaxis left
                else
                    yyaxis right
                end
                plot(obj.sweepValues,obj.solution(compIndex(n),:),lineSpec);
            end
            
            subplot(1,1,1);
            title('DC Sweep');
            legend(obj.solvedValuesNames(compIndex));
            grid on;
            xlabel(obj.componentList(obj.inputComponentCompNum).getUnit);
            yyaxis left
            ylabel('Voltage (V)');
            yyaxis right
            ylabel('Current (A)');
        end
        
        function [f] = outputAC(obj,compIndex)
            f = figure;
            for n = 1:length(compIndex)
                lineSpec = circ.getLineSpec(n);
                
                G = obj.solution(compIndex(n),:)./obj.acSweepInput;
                subplot(2,2,1)
                hold on;
                plot(obj.sweepValues,20*log10(abs(G)),lineSpec);
                
                subplot(2,2,2)
                hold on;
                plot(obj.sweepValues,abs(obj.solution(compIndex(n),:)),lineSpec);
                
                subplot(2,2,3)
                hold on;
                plot(obj.sweepValues,rad2deg(angle(G)),lineSpec);
                
                subplot(2,2,4)
                hold on;
                plot(obj.sweepValues,rad2deg(angle(obj.solution(compIndex(n),:))),lineSpec);
            end
            subplot(2,2,1)
            xlabel('Frequency (Hz)');
            ylabel('Gain (dB)');
            title('Bode Plot');
            legend(obj.solvedValuesNames(compIndex));
            grid on
            
            subplot(2,2,2)
            xlabel('Frequency (Hz)');
            ylabel('Output Voltage (V)');
            title('AC Sweep');
            legend(obj.solvedValuesNames(compIndex));
            grid on
            
            subplot(2,2,3)
            xlabel('Frequency (Hz)');
            ylabel('Phase (deg)');
            grid on
            
            subplot(2,2,4)
            xlabel('Frequency (Hz)');
            ylabel('Phase (deg)');
            grid on
        end
        
        function [h] = outputTransient(obj,compIndex)
            h = figure;
            for n = 1:length(compIndex)
                lineSpec = circ.getLineSpec(n);
                
                subplot(2,1,1);
                hold on;
                if obj.isVarVoltage(compIndex(n))
                    yyaxis left
                else
                    yyaxis right
                end
                plot(obj.sweepValues,obj.solution(compIndex(n),:),lineSpec);
                
                subplot(2,1,2);
                hold on;
                if obj.isVarVoltage(compIndex(n))
                    yyaxis left
                else
                    yyaxis right
                end
                [f,S] = circ.timeToFreq(obj.sweepValues,obj.solution(compIndex(n),:));
                plot(f,S,lineSpec);
                
            end
            subplot(2,1,1);
            title('Transient Sweep');
            legend(obj.solvedValuesNames(compIndex));
            grid on;
            xlabel('Time (s)');
            yyaxis left
            ylabel('Voltage (V)');
            yyaxis right
            ylabel('Current (A)');
            
            subplot(2,1,2);
            title('Frequency Spectrum');
            legend(obj.solvedValuesNames(compIndex));
            grid on;
            xlabel('Frequency (Hz)');
            yyaxis left
            ylabel('Voltage (V)');
            yyaxis right
            ylabel('Current (A)');
            xlim([0,50])
        end
        
        function [compIndex] = compToIndex(obj,comps)
            if length(comps) == 1 && isequal(comps{1},'all')
                compIndex = 1:obj.numEqn;
                return;
            end
            compNum = length(comps);
            compIndex = zeros(1,compNum);
            for n = 1:compNum
                compIndex(n) = obj.getCompNum(comps{n});
                if compIndex(n) == -1
                    fprintf('%s is an invalid component\n',comps{n});
                end
            end
            compIndex(compIndex == -1) = [];
        end
        
        %%General
        
        function [c] = getComp(obj,name)
            c = obj.componentList(obj.getCompNum(name));
        end
        
        function [obj] = setupIfNotAlready(obj)
            if obj.status == 0
                obj = obj.setup();
            end
        end
        
        function [tf] = isVarVoltage(obj,index)
            currName = obj.solvedValuesNames{index};
            tf = currName(1) == 'V';
        end
        
        function [out] = varUnit(obj,index)
            if obj.isVarVoltage(index)
                out = 'V';
            else
                out = 'A';
            end
        end
        
        function [out] = getCompNum(obj,name)
            if ~ischar(name)
                out = name + obj.numComp;
                if out > obj.numEqn
                    out = -1;
                end
                return;
            end
            
            for currComp = obj.componentList
                if (length(currComp.name) == length(name)) & (currComp.name == name)
                    out = currComp.num;
                    return;
                end
            end
            out = -1;
        end
    end
    
    methods
        %%Initialize
        function obj = circ()
            obj.componentList = [];
            obj.numComp = 0;
            obj.voltageNodes = [];
            obj.status = 0;
        end
        
        %%Input
        
        function obj = add(obj,name,n1,n2,val,rel)
            if sum(obj.voltageNodes == n1) == 0
                obj.voltageNodes = [obj.voltageNodes,n1];
            end
            if sum(obj.voltageNodes == n2) == 0
                obj.voltageNodes = [obj.voltageNodes,n2];
            end
            obj.numComp = obj.numComp + 1;
            if nargin == 5
                obj.componentList = [obj.componentList,comp(name,n1,n2,val,obj.numComp)];
            else
                obj.componentList = [obj.componentList,comp(name,n1,n2,val,obj.numComp,rel)];
            end
            obj.status = 0;
        end
        
        function obj = change(obj,name,newVal)
            compNum = obj.getCompNum(name);
            obj.componentList(compNum).val = newVal;
            obj = assignVal(obj,obj.componentList(compNum));
        end
        
        function [obj] = setIn(obj,name)
            obj.inputComponentName = name;
        end
        
        %%Setup
        
        function [obj] = setup(obj)
            obj.voltageNodes = sort(obj.voltageNodes);
            obj.numNodes = length(obj.voltageNodes)-1;
            obj.numEqn = obj.numNodes + obj.numComp;
            obj.solvedValuesNames = cell(obj.numEqn,1);
            obj.Gmat = zeros(obj.numEqn,obj.numEqn);
            obj.Cmat = zeros(obj.numEqn,obj.numEqn);
            obj.Fmat = zeros(obj.numEqn,1);
            obj = obj.performKCL();
            obj = obj.analyzeComps();
            obj.inputComponentCompNum = obj.getCompNum(obj.inputComponentName);
            obj.inputComponentEqnNum = obj.inputComponentCompNum + obj.numNodes;
            obj.status = 1;
        end
        
        %%Solve
        
        function [obj] = opPoint(obj,Vin,f)
            obj = obj.setSweepValues(Vin);
            if nargin == 2
                obj = obj.solve(Vin);
            else
                obj.sweepValues = [obj.sweepValues,f];
                obj = obj.solve(Vin,f);
            end
            obj.solution = obj.Vmat;
            obj.status = 2;
        end
        
        function [obj] = DCsweep(obj,sweepValues)
            obj = obj.setSweepValues(sweepValues);
            for n = 1:obj.numSweepValues
                obj = obj.solve(sweepValues(n));
                obj.solution(:,n) = obj.Vmat;
            end
            obj.status = 3;
        end
        
        function [obj] = ACsweep(obj,sweepValues,acSweepInput)
            if nargin == 2
                acSweepInput = 1;
            end
            obj.acSweepInput = acSweepInput;
            obj = obj.setSweepValues(sweepValues);
            for n = 1:obj.numSweepValues
                obj = obj.solve(acSweepInput,2*pi*sweepValues(n));
                obj.solution(:,n) = obj.Vmat;
            end
            obj.status = 4;
        end
        
        function [obj] = transient(obj,sweepValues,inputV,randElName,randElMean,randElStd)
            if nargin == 3
                varyEl = false;
            else
                varyEl = true;
            end
            
            obj = obj.setSweepValues(sweepValues);
            tmp = obj.Fmat;
            obj.Fmat = repmat(obj.Fmat, [1,obj.numSweepValues]);
            obj.Fmat(obj.inputComponentEqnNum,:) = inputV;
            if varyEl
                obj.Fmat(obj.getCompNum(randElName),:) = normrnd(randElMean,randElStd,1,obj.numSweepValues);
            end
            obj.solution(:,1) = obj.Fmat(:,1);
            for n = 2:obj.numSweepValues
                dt = obj.sweepValues(n) - obj.sweepValues(n-1);
                obj.solution(:,n) = mldivide(obj.Cmat./dt + obj.Gmat,obj.Cmat*obj.solution(:,n-1)./dt + obj.Fmat(:,n-1));
            end
            obj.Fmat = tmp;
            obj.status = 5;
        end
        
        
        %%Output
        
        function out = I(obj,num)
            if ischar(num)
                out = obj.solution(obj.getCompNum(num),:);
            else
                out = obj.solution(num,:);
            end
        end
        
        function out = V(obj,num)
            out = obj.solution(obj.getCompNum(num),:);
        end
        
        function [] = printStatus(obj)
            switch obj
                case 0
                    fprintf('Not Initialized\n');
                case 1
                    fprintf('Setuped\n');
                case 2
                    fprintf('Solved at a DC operating Point\n');
                case 3
                    fprintf('DC Sweep Performed\n');
                case 4
                    fprintf('AC Sweep Performed\n');
                case 5
                    fprintf('Transiant Analysis Performed\n');
            end
        end
        
        function [f] = output(obj,varargin)
            compIndex = compToIndex(obj,varargin);
            f = -1;
            switch obj.status
                case 2
                    obj.outputOpPoint(compIndex);
                case 3
                    f = obj.outputDC(compIndex);
                case 4
                    f = obj.outputAC(compIndex);
                case 5
                    f = obj.outputTransient(compIndex);
                otherwise
                    return;
            end
        end
        
        function [out] = getSol(obj,varargin)
            compIndex = compToIndex(obj,varargin);
            out = obj.solution(compIndex,:);
        end
        
        function [] = printSol(obj)
            if obj.status < 2
                return;
            end
            for n = 1:obj.numEqn
                if imag(obj.Vmat(n)) == 0
                    fprintf('%s = %f\n',obj.solvedValuesNames{n},obj.Vmat(n));
                else
                    fprintf('%s = %f+%fj\n',obj.solvedValuesNames{n},real(obj.Vmat(n)),imag(obj.Vmat(n)));
                end
            end
        end
        
        function [] = printEqn(obj)
            for n = 1:obj.numEqn
                ls = '';
                rs = '';
                text = sprintf('%.1f',abs(obj.Fmat(n)));
                [ls,rs] = circ.addToSide(ls,rs,text,-obj.Fmat(n));
                for m = 1:obj.numEqn
                    if obj.Gmat(n,m) ~= 0
                        text = sprintf('%s%s',circ.printVal(obj.Gmat(n,m)),obj.solvedValuesNames{m});
                        [ls,rs] = circ.addToSide(ls,rs,text,obj.Gmat(n,m));
                    end
                    if obj.Cmat(n,m) ~= 0
                        text = sprintf('%sjw%s',circ.printVal(obj.Cmat(n,m)),obj.solvedValuesNames{m});
                        [ls,rs] = circ.addToSide(ls,rs,text,obj.Cmat(n,m));
                    end
                end
                if isempty(ls)
                    ls = '0';
                end
                if isempty(rs)
                    rs = '0';
                end
                fprintf('%s = %s\n',ls,rs);
            end
            fprintf('\n');
        end
        
        function [] = printMat(obj,matRep,latex)
            if nargin == 2
                latex = false;
            end
            obj = setupIfNotAlready(obj);
            switch lower(matRep(1))
                case 'g'
                    mat = obj.Gmat;
                case 'c'
                    mat = obj.Cmat;
                case 'f'
                    mat = obj.Fmat;
            end
            fprintf(circ.makeLength('Eqn',8));
            
            if size(mat,2) == 1
                if latex
                        fprintf('& \\\\');
                 end
            else
                for n = 1:obj.numEqn
                    if latex
                        fprintf('& ');
                    end
                    fprintf('%s',circ.makeLength(sprintf(obj.solvedValuesNames{n}),8));
                end
                if latex
                    fprintf(' \\\\');
                end
            end
            fprintf('\n');
            for n = 1:obj.numEqn
                fprintf('%s',circ.makeLength(sprintf('%i',n),8));
                for m = 1:size(mat,2)
                    if latex
                        fprintf('& ');
                    end
                    fprintf('%s',circ.formatVal(mat(n,m),8));
                end
                if latex
                    fprintf(' \\\\');
                end
                fprintf('\n');
            end
            fprintf('\n');
        end
        
    end
end

