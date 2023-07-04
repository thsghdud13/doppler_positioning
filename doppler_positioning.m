%% gps 응용 23 도플러 변이 값으로 수신기 위치 계산
clc; close all; clearvars
eph = ReadEPH('brdc3530.22n');
[Lat, Lon, TEC] = ReadGIM('igsg3530.22i');
QMfile = 'QM_hangangstop221219_1816';
arrQM = load(QMfile);
TruePos_LLH = [37.522531765 126.940182578 16.082];
TruePos_XYZ = gd2xyz(TruePos_LLH);
gw = 2241;
%% 상수, 변수 정의
CCC = 299792458; % 빛의 속도
obsType = 103;
lmbd1 = CCC/(1575.42*10^6); % L1 주파수의 파장
%% QM 선별
QM_pre = arrQM(arrQM(:,3) == obsType,:);
QM = QM_pre(QM_pre(:,2) < 199,:);
FinalTTs = unique(QM(:,1));
%% 추정에 필요한 반복 조건 및 초기값 설정
MaxIter = 16;
x = [TruePos_XYZ 0]; 
x = x';
%% 추정과정 시작
NoEpochs = length(FinalTTs);
estm = zeros(NoEpochs, 5); % 결과 저장 행렬
nEst = 0; % 추정이 안될 경우 -> 인덱스를 증가시키기 위해서
NEV = zeros(NoEpochs,3);
Error_3d = zeros(NoEpochs,1);
Error_clock = zeros(NoEpochs,1);
for kE = 1:NoEpochs
    indx1e = find(QM(:,1) == FinalTTs(kE));
    QM1e = QM(indx1e,:);
    gs = QM1e(1,1);
    [r,c] = size(QM1e);
    NoSats = r;
    if(NoSats < 5)
        continue;
    end
    snr=QM1e(:,7);
    vec_site = x(1:4);
    dop_epoch = QM1e(:,6);
    dop_epoch = dop_epoch.*lmbd1;
    B = zeros(NoSats,1);
    B2 = zeros(4,1);
    A2 = zeros(NoSats,NoSats);
    Amat = zeros(NoSats,4);
    elaz=zeros(NoSats,2);
    randomError = rand(1)*1000;
    for kI = 1:MaxIter
        vec_pos = vec_site(1:3)';
        w=zeros(NoSats,NoSats);
        ww=zeros(NoSats,NoSats);
        ZHD = TropGPTh(vec_pos,gw,gs);
        for kS = 1:NoSats
            PRN = QM1e(kS, 2);
            if PRN > 199
                continue;
            end

            %% 신호세기 가중치

            for sats = 1 : NoSats
                if snr(sats,1)<30
                    ww(sats,sats)=0;
                else
                    ww(sats,sats)=snr(sats,1)/10;
                end
            end



            obs = QM1e(kS,4);
            STT = obs/CCC;
            tCorr = gs - STT;
            ieph = PickEPH(eph,PRN,tCorr);
            vec_sat = getSatPos(eph(ieph,:),tCorr);
            vec_sat = RotSatPos(vec_sat,STT);
            vec_sat_before = getSatPos(eph(ieph,:),tCorr-0.1);
            vec_sat_before = RotSatPos(vec_sat_before, STT);
            vec_sat_after = getSatPos(eph(ieph,:),tCorr+0.1);
            vec_sat_after = RotSatPos(vec_sat_after, STT);
            velo_sat = getSat_D(eph,ieph,gs);
            velo_sat_after = getSat_D(eph,ieph,gs+0.1);
            velo_sat_before = getSat_D(eph,ieph,gs-0.1);

            
            elaz(kS,:)=getElAz(x(1:3)',vec_sat); % 고도각 가중치
            w(kS,kS)=2*sin(deg2rad((elaz(kS,1))));


            Amat(kS,1:3) = velo_sat; % A행렬
            Amat(kS,4) = norm(vec_sat);
            Y = vec_site(1:3) - vec_sat'; % 위성에서 수신기로의 벡터
            Ym = sqrt(Y'*Y);
            %% 오차 모델링
            eph_sat = eph(ieph,:);
            b = eph_sat(4);
            sat_drift = b; % 위성시계 오차 변화율
            dTrop = dTropDot(vec_sat_before,vec_sat,vec_sat_after,gw,tCorr,vec_pos,ZHD); % 대류권오차 변화율
            dIono = IONO(Lat,Lon,TEC,vec_sat_before,vec_sat,vec_sat_after,vec_pos,tCorr); % 이온층오차 변화율
            dSagnac = dRot(velo_sat,vec_site); % 지구자전효과 변화율
            dRel  = relativistic(vec_sat,vec_sat_after,velo_sat,velo_sat_after); % 상대성효과 변화율
            dop_tmp = dop_epoch(kS); % 관측된 의사거리속도
            B(kS) = vec_sat * velo_sat' + (dop_tmp)*Ym + (vec_site(4)) * (Amat(kS,4) - Ym);% + (Ym)*(CCC*sat_drift - dTrop  - dSagnac - dIono - CCC*dRel);
        end
        B2 = Amat'*ww* B;
        A2 = Amat'*ww* Amat;
        vec_site = A2\B2;
    end

    nowPos = vec_site(1:3);
    nowPos = nowPos';
    nowGd = xyz2gd(nowPos);
    fprintf("\n\n");
    dXYZ = nowPos - TruePos_XYZ;
    fprintf("3dError : %10.10f\n",norm(dXYZ));
    dNEV = xyz2topo(dXYZ,TruePos_LLH(1),TruePos_LLH(2));
    NEV(kE,:) = dNEV;
    Error_3d(kE,:) = norm(dXYZ);
    fprintf("dN: %10.10f dE: %10.10f dV: %10.10f\n",dNEV(1),dNEV(2),dNEV(3));
    fprintf("plane Error : %10.10f\n", sqrt(dNEV(1).^2 + dNEV(2).^2));
    estm(kE,1) = gs;
    estm(kE,2:4) = nowPos;
    estm(kE,5) = vec_site(4);
    Error_clock(kE,:) = vec_site(4);
    aaa = gs2h24(gs);
end

dN = NEV(:,1);
dE = NEV(:,2);
dV = NEV(:,3);
Nrms = rms(dN);
Erms = rms(dE);
Vrms = rms(dV);
Hor = norm([Nrms,Erms]);
rms_3d = norm([Nrms, Erms, Vrms]);

fprintf("Hor Rms = %10.10f \n",Hor);
fprintf("3d Rms = %10.10f \n",rms_3d);
fprintf("V Rms = %10.10f \n",Vrms);

% % subplot(3,4,[5,6]);
% plot(dE,dN,'ob');
% axis([-3000 3000 -3000 3000]);
% % axis tight;
% grid on; ylabel('dN[m]'); xlabel('dE[m]');

% subplot(3,4,[1,2]);
% plot(gs2h24(estm(:,1)),dN,'om');
% axis tight;
% grid on; ylabel('dN[m]'); xlabel('hours');
% 
% 
% subplot(3,4,[9,10]);
% plot(gs2h24(estm(:,1)),dE,'om');
% axis tight;
% grid on; ylabel('dE[m]'); xlabel('hours');
% 
% subplot(3,4,[3,4])
plot(gs2h24(estm(:,1)),Error_3d,"or");
% hold on;
% plot(gs2h24(estm(:,1)),NumSat*1000,"og")
% xlim([0 24])
% grid on;
% ylabel('3D Error[m]'); xlabel('hours');
% 
% subplot(3,4,[7,8])
% plot(gs2h24(estm(:,1)),Error_clock,"or");
% % xlim([0 24])
% grid on;
% ylabel('delta w'); xlabel('hours');
% 
% 
% % subplot(3,4,[11 12])
% plot(gs2h24(estm(:,1)),dV,"or");
% % xlim([0 24])
% grid on;
% ylabel('dV[m]'); xlabel('hours');





