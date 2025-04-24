/*
 * SignalBrowser.c
 *
 *  Created on: 7 мар. 2024 г.
 *      Author: ES
 */

#include "SignalBrowser.h"
#include "pelengator.h"
#include "modulations.h"
 #include "axi_adc_core.h"
 #include "axi_dmac.h"
 #include "app_talise.h"
 #include "talise_error_types.h"
 #include "utils.h"
#include "xtime_l.h"
#include "registers.h"
#include "utils.h"
#include "crc32.h"
#include <stdio.h>
#include "xil_printf.h"

#include "math.h"

// Количество отчетов
 // 32 бита = 4 байта

//// Буфер для хранения отчетов
//#define FILTER_ORDER 6
//#define NUM_SECTIONS (FILTER_ORDER / 2)

volatile u8 buf1[64];

u16 local_count_buff=0;

u16 local_count_buff2=0;


///Коефициенты под фильтр и децимацию
#define N1 65536  // количество отсчетов сигнала 8192/16384//32768//65536//131072
#define FC 30e6 // частота среза фильтра, 30 МГц
#define FS 245e6 // частота дискретизации, 245 МГц


volatile comp fftpsd_rx1[256];
volatile comp fftpsd_rx2[256];
volatile float  BUF_ADC_I[N1];
volatile float  BUF_ADC_Q[N1];
volatile s16 BUF_ADC_Q_filter[N1];
volatile s16 BUF_ADC_I_filter[N1];

volatile double phase_unwrapped_deg_1[N1];
volatile double jump_indices_1[N1];


volatile u32 adcpsd_rx1[256];
volatile u32 adcpsd_rx2[256];

volatile u16 global_ind = 0;
volatile u8 global_AntCount=0;
volatile u8 global_FreqRepeatCount=0;
volatile u8 global_MaxHoldRepeatCount=0;

volatile double Pob_min_array[128];
void init_HardwareTask()
{


	Pelengator.aiFmin[0] = 500;//2490;//2250;
	Pelengator.aiFCount[0] = 74;//74;//23
	Pelengator.AntNum = 11;
	Pelengator.RxGain_Array[0] = 255;
	Pelengator.MaxThreshold = -50;

	Pelengator.aDelayCycle = 30;
	Pelengator.MaxHoldDelay = 80;
	Pelengator.MaxHoldRepeatCount =20;

/*	Pelengator.aDelayCycle = 0;
	Pelengator.MaxHoldDelay = 4000;
	Pelengator.MaxHoldRepeatCount =1;
*/


	Pelengator.sensivity = 5.0;
	Pelengator.max_delta_samples = 3;//3
	Pelengator.max_delta_rf = 10.0;//5.0;

	Pelengator.porog_fz = 1;
	Pelengator.porog_band = 0.6;


	Pelengator.Len=256;
	Pelengator.StartAnalyse=3;
	Pelengator.StopAnalyse=253;

	Pelengator.TensCeiling=20;
	Pelengator.Absence=20; //
	Pelengator.AbsenceTime= 20;

	Pelengator.DirectRXcA = 0xf;
	Pelengator.DirectRXcD0 = 1;//0 enable test gen
	Pelengator.DirectRXcD1 = 1;
	Pelengator.DirectRXcATT = 0x7;
	Pelengator.DirectRXcDD = 4;
	Pelengator.DirectRXcAA = 0;

	volatile u8 iA = 0;
	volatile u16 i=0;
	for(i=15;i>15-11;i--)
	{
		Pelengator.RXcA[iA++] = i;
	}
	//for(i=0;i<11;i++) xil_printf("Pelengator.RXcA[%d]=%d \r\n",i, Pelengator.RXcA[i]);
	for(i=0;i<128;i++)
	{
		Pelengator.RXcATT[i] = 0x7;
		Pelengator.RXcDD[i] = 4;
		Pelengator.RXcAA[i] = 0;
	}

	u16 conn_index=0;
	u64 tmp = Pelengator.aiFmin[0];
	//xil_printf("start tmp = %d \r\n",(u32) tmp);
    for(i = 0;i<Pelengator.aiFCount[0];i++)
    {
        //u64 tmp = Pelengator.aiFCount[i];

        Pelengator.RXcATT[conn_index]=0x7;
        if ( (tmp > 0) && (tmp <= 6e3) )
        {
            Pelengator.RXcDD[conn_index] = 0x4;
			Pelengator.RXcAA[conn_index] = 0x0;
        }
        else if ( (tmp > 6e3) && (tmp < 10e3) )
        {
            Pelengator.RXcDD[conn_index] = 0x3;
            Pelengator.RXcAA[conn_index] = 0x1;
        }
        else if ( (tmp >= 10e3) && (tmp < 14e3) )
        {
            Pelengator.RXcDD[conn_index] = 0x1;
            Pelengator.RXcAA[conn_index] = 0x2;
        }
        else if ( (tmp >= 14e3) && (tmp < 20e3) )
        {
            Pelengator.RXcDD[conn_index] = 0x2;
            Pelengator.RXcAA[conn_index] = 0x4;
            //xil_printf("RXcDD=%d RXcAA=%d tmp = %d\r\n", Pelengator.RXcDD[conn_index], Pelengator.RXcAA[conn_index], (u32)tmp);
        }
        tmp += (Pelengator.CalibrationFreqStep/1e6);


        conn_index++;
    }

   // xil_printf("end tmp = %d conn_index =%d\r\n",(u32)tmp, conn_index);


	xil_printf("start freq = %d end freq = %d \r\n",(u32) (Pelengator.aiFmin[0]-Pelengator.CalibrationFreqStep/2e6), (u32) (Pelengator.aiFmin[0]+Pelengator.aiFCount[0]*Pelengator.CalibrationFreqStep/1e6-Pelengator.CalibrationFreqStep/2e6));
	u32 DataSize = Create_HardwareTask(buf1);
	xil_printf("Task DataSize=%d byte\r\n", DataSize);
	Xil_DCacheInvalidateRange((UINTPTR)(RX_BUFFER_BASE),  DataSize);

	//TSignal S;
	//memset(&S, 0, sizeof(TSignal));

	//xil_printf("sizeof(fftpsd_rx1)=%d\r\n", sizeof(fftpsd_rx1));

	//memset(&fftpsd_rx1, 0, sizeof(fftpsd_rx1));
	//memset(&fftpsd_rx2, 0, sizeof(fftpsd_rx2));

	/*volatile u8 i=0;

	for(i=0;i<40;i++)
	{
		S.Freq = i*10;
		S.fc_max = (u16)round(S.Freq);
		S.power=-(i*2);
		AddToSignalTable(S);
	}
	*/
/*
	S.Freq = 1;
	S.fc_max = (u16)round(S.Freq);
	S.power=-52;
	AddToSignalTable(S);

	S.Freq = 20;
	S.fc_max = (u16)round(S.Freq);
	S.power=-50;
	AddToSignalTable(S);

	S.Freq = 40;
	S.fc_max = (u16)round(S.Freq);
	S.power=-54;
	AddToSignalTable(S);

	SortByDecrease(0);
	print_SignalTable();
*/


	return;
}
void main_HardwareVSKTask()
{
	xil_printf("void main_HardwareVSKTask()\r\n");


	Pelengator.aiFmin[0] = 4010;//2490;//2250;
	Pelengator.aiFCount[0] = 1;
	Pelengator.aDelayCycle = 10;
	Pelengator.MaxHoldDelay = 4000;
	Pelengator.MaxHoldRepeatCount = 1;
	Pelengator.Freq = 4000;
	Pelengator.AntNum = 11;
	Pelengator.RxGain_Array[0] = 255;
	Pelengator.MaxThreshold = -80;

	Pelengator.sensivity = 3.0;
	Pelengator.max_delta_samples = 2;//3
	Pelengator.max_delta_rf = 50.0;



	Pelengator.Len=256;
	Pelengator.StartAnalyse=3;
	Pelengator.StopAnalyse=253;

	Pelengator.TensCeiling=1;
	Pelengator.Absence=1; //
	Pelengator.AbsenceTime= 10;




	Pelengator.DirectRXcA = 0xf;
	Pelengator.DirectRXcD0 = 0;//0 enable test gen
	Pelengator.DirectRXcD1 = 0;
	Pelengator.DirectRXcATT = 0x7;
	Pelengator.DirectRXcDD = 4;
	Pelengator.DirectRXcAA = 0;

	volatile u8 iA = 0;
	volatile u8 i=0;
	for(i=15;i>15-11;i--)
		Pelengator.RXcA[iA++] = i;


	for(u8 i=0;i<128;i++)
	{
		Pelengator.RXcATT[i] = 0x7;
		Pelengator.RXcDD[i] = 4;//3;//4;
		Pelengator.RXcAA[i] = 0;//1;//0;
	}





	//xil_printf("start tmp = %d \r\n",(u32) tmp);
   /* for(i = 0;i<Pelengator.aiFCount[0];i++)
    {
        //u64 tmp = Pelengator.aiFCount[i];

        Pelengator.RXcATT[conn_index]=0x7;
        if ( (tmp > 0) && (tmp <= 6e3) )
        {
            Pelengator.RXcDD[conn_index] = 0x4;
			Pelengator.RXcAA[conn_index] = 0x0;
        }
        else if ( (tmp > 6e3) && (tmp < 10e3) )
        {
            Pelengator.RXcDD[conn_index] = 0x3;
            Pelengator.RXcAA[conn_index] = 0x1;
        }
        else if ( (tmp >= 10e3) && (tmp < 14e3) )
        {
            Pelengator.RXcDD[conn_index] = 0x1;
            Pelengator.RXcAA[conn_index] = 0x2;
        }
        else if ( (tmp >= 14e3) && (tmp < 20e3) )
        {
            Pelengator.RXcDD[conn_index] = 0x2;
            Pelengator.RXcAA[conn_index] = 0x4;
            //xil_printf("RXcDD=%d RXcAA=%d tmp = %d\r\n", Pelengator.RXcDD[conn_index], Pelengator.RXcAA[conn_index], (u32)tmp);
        }
        tmp += (Pelengator.CalibrationFreqStep/1e6);


        conn_index++;
    }*/



	SetDirectExternalPorts();
	u32 DataSize = Create_HardwareTask(buf1);





	Xil_DCacheInvalidateRange((UINTPTR)(RX_BUFFER_BASE),  DataSize);




	Cmd11_type2xScanWBTask_hardware_scan(buf1);
	FindMax_Scan(buf1);

	Pelengator.VSK_rezult = 0;
	double vsk_freq = 4000;
	//if(Pelengator.SignalCnt == 0)
	xil_printf("Pelengator.SignalCnt() = %d \r\n", Pelengator.SignalCnt);
	if(Pelengator.SignalCnt==0)Pelengator.VSK_rezult = ADSINGLE_VSK_STATUS_FAIL;
	for(u8 i=0;i<Pelengator.SignalCnt;i++)
	{
  	  TSignal * TmpSignal=&(SignalTable[i]);//(sSignal*)SignalList->Items[i];

  	  float Tmp = fabs(vsk_freq - TmpSignal->Freq);
  	//xil_printf("TmpSignal[%d]->Freq=%d\r\n",i,(s32)TmpSignal->Freq);
  	  if (Tmp < 2)
  	{
  		//xil_printf("TmpSignal[%d]->Freq=%d\r\n",i,(s32)TmpSignal->Freq);
  		//xil_printf("Tmp=%d\r\n",(s32)Tmp);
  		  Pelengator.VSK_rezult = ADSINGLE_VSK_STATUS_SUCCESS;

  		  break;
  	}
  	  else Pelengator.VSK_rezult = ADSINGLE_VSK_STATUS_FAIL;
	}
	if(Pelengator.Start_VSK_rezult == ADSINGLE_VSK_STATUS_NO)
		Pelengator.Start_VSK_rezult = Pelengator.VSK_rezult;
	CalculatePeleng();
	//print_SignalTable();
xil_printf("vsk end\r\n");

}
void main_HardwareTask()
{
	xil_printf("main_HardwareTask \r\n");


	Cmd11_type2xScanWBTask_hardware_scan(buf1);
	FindMax_Scan(buf1);
	//xil_printf("CheckSignalsTime \r\n");
	CheckSignalsTime();
	//xil_printf("end task \r\n");
	//SortByDecrease(0);
	// CalculatePeleng();
	//print_SignalTable();


}
u32 Create_HardwareTask(u8* buf)
{
	//xil_printf("Create_HardwareTask \r\n");
	volatile u8 *DataBuf = (buf + 20); // (buf + 20);
	*(u16*)(buf + PAYLOAD_CMD_OFFSET) = 0; // CmdNum

	(*(u16 *)(DataBuf + 0)) = 512; //splencode

	DataBuf[4] = 1; // aFDcount
	volatile u8 aFDcount = DataBuf[4];

	(*(u16 *)(DataBuf + 5)) = 0; // aTimerPeriod

	(*(u16 *)(DataBuf + 7)) = Pelengator.aDelayCycle;// 4; // aDelayCycle

	volatile u16 aDetectThreshold = (*(u16 *)(DataBuf + 9));
	(*(u16 *)(DataBuf + 9)) = 0;// aDetectThreshold
	(*(u8 *)(DataBuf + 11 )) = Pelengator.AntNum; // aiAntCount;
	(*(u8 *)(DataBuf + 12)) = 1;// aiTaskNumber
	(*(u8 *)(DataBuf + 13)) =0 ; // aiExternalSyncEn

	volatile u16 i = 0;
	volatile u16 offset = 14;
	for(i = 0;i<aFDcount;i++)
	{
		(*(u16 *)(DataBuf + offset+0 + i*5)) = Pelengator.aiFmin[i];
		(*(u16 *)(DataBuf + offset+2 + i*5)) = Pelengator.aiFCount[i];
		(*(u8  *)(DataBuf + offset+4 + i*5)) = Pelengator.RxGain_Array[i];
	}

	   volatile u8 TaskRepeatCount;// = (*(u8 *)(DataBuf + offset + aFDcount*(5)));
	   (*(u8 *)(DataBuf + offset + aFDcount*(5))) = 1; // TaskRepeatCount
	   //xil_printf("TaskRepeatCount = %d \r\n", TaskRepeatCount);
	   volatile u8 FreqRepeatCount;// = (*(u8 *)(DataBuf + offset + aFDcount*(5) + 1));
	   (*(u8 *)(DataBuf + offset + aFDcount*(5) + 1)) = 1; //FreqRepeatCount
	   (*(u8 *)(DataBuf + offset + aFDcount*(5) + 2)) = Pelengator.MaxHoldRepeatCount; // 1;//MaxHoldRepeatCount
	   (*(u16 *)(DataBuf + offset + aFDcount*(5) + 3)) = Pelengator.MaxHoldDelay;// MaxHoldDelay;
	   (*(u16 *)(DataBuf + offset + aFDcount*(5) + 5)) = 0;// DetectFreq;
	   (*(u16 *)(DataBuf + offset + aFDcount*(5) + 7)) = 0;// DetectFreqDelta;

	   u32 DataSize = 0;

	   //for(FDcountINDEX=0;FDcountINDEX < aFDcount; FDcountINDEX++)
	   {
		   DataSize = ((Pelengator.Len << Pelengator.ch_num_shift)) * Pelengator.aiFCount[0] * Pelengator.AntNum * Pelengator.MaxHoldRepeatCount; // FreqRepeatCount  * MaxHoldRepeatCount * aiAntCount;
	   }

	return DataSize;
}
u32 Create_HardwareVSKTask(u8* buf)
{
	//xil_printf("Create_HardwareTask \r\n");
	volatile u8 *DataBuf = (buf + 20); // (buf + 20);
	*(u16*)(buf + PAYLOAD_CMD_OFFSET) = 0; // CmdNum

	(*(u16 *)(DataBuf + 0)) = 512; //splencode

	DataBuf[4] = 1; // aFDcount
	volatile u8 aFDcount = DataBuf[4];

	(*(u16 *)(DataBuf + 5)) = 0; // aTimerPeriod

	(*(u16 *)(DataBuf + 7)) = Pelengator.vsk_aDelayCycle;// 4; // aDelayCycle

	volatile u16 aDetectThreshold = (*(u16 *)(DataBuf + 9));
	(*(u16 *)(DataBuf + 9)) = 0;// aDetectThreshold
	(*(u8 *)(DataBuf + 11 )) = Pelengator.AntNum;//11; // aiAntCount;
	(*(u8 *)(DataBuf + 12)) = 1;// aiTaskNumber
	(*(u8 *)(DataBuf + 13)) =0 ; // aiExternalSyncEn

	volatile u16 i = 0;
	volatile u16 offset = 14;
	//for(i = 0;i<aFDcount;i++)
	{
		(*(u16 *)(DataBuf + offset+0 + i*5)) = Pelengator.vsk_aiFmin;
		(*(u16 *)(DataBuf + offset+2 + i*5)) = Pelengator.vsk_aiFCount;
		(*(u8  *)(DataBuf + offset+4 + i*5)) = Pelengator.RxGain_Array[0];
	}

	   volatile u8 TaskRepeatCount;// = (*(u8 *)(DataBuf + offset + aFDcount*(5)));
	   (*(u8 *)(DataBuf + offset + aFDcount*(5))) = 1; // TaskRepeatCount
	   //xil_printf("TaskRepeatCount = %d \r\n", TaskRepeatCount);
	   volatile u8 FreqRepeatCount;// = (*(u8 *)(DataBuf + offset + aFDcount*(5) + 1));
	   (*(u8 *)(DataBuf + offset + aFDcount*(5) + 1)) = 1; //FreqRepeatCount
	   (*(u8 *)(DataBuf + offset + aFDcount*(5) + 2)) = Pelengator.MaxHoldRepeatCount;//1;//MaxHoldRepeatCount
	   (*(u16 *)(DataBuf + offset + aFDcount*(5) + 3)) =Pelengator.MaxHoldDelay;// MaxHoldDelay;
	   (*(u16 *)(DataBuf + offset + aFDcount*(5) + 5)) = 0;// DetectFreq;
	   (*(u16 *)(DataBuf + offset + aFDcount*(5) + 7)) = 0;// DetectFreqDelta;

	   u32 DataSize = 0;

	   //for(FDcountINDEX=0;FDcountINDEX < aFDcount; FDcountINDEX++)
	   {
		   DataSize = ((Pelengator.Len << Pelengator.ch_num_shift)) * Pelengator.aiFCount[0] * 11; // FreqRepeatCount  * MaxHoldRepeatCount * aiAntCount;
	   }

	return DataSize;
}
u32 Cmd11_type2xScanWBTask_hardware_scan(u8* buf)
{
	xil_printf("Cmd11_type2xScanWBTask_hardware_scan \r\n");

	//XGpioPs_WritePin(&Gpio, 10, 0x1);
	//XGpioPs_WritePin(&Gpio, 12, 0x1);
	//DiscreteSet(rCPU_LEDs, 1);

	//RegDiscreteSet(rLEDS,2);


	volatile u16 CmdNum = *(u16*)(buf + PAYLOAD_CMD_OFFSET);
	//xil_printf("Cmd11_type2xScanWBTask Cmd %d rcv!\n\r",CmdNum);


	volatile u8 *DataBuf = (buf + 20);
	volatile u32 splencode = (*(u16 *)(DataBuf + 0));
	splencode &= ~(0x3);
	splencode >>= 1;

	volatile u8 aFDcount = DataBuf[4];

	volatile u16 aTimerPeriod = (*(u16 *)(DataBuf + 5));
	volatile u16 aDelayCycle = (*(u16 *)(DataBuf + 7));

	volatile u16 aDetectThreshold = (*(u16 *)(DataBuf + 9));
	Pelengator.Threshold = aDetectThreshold;

	volatile u8 aiAntCount = (*(u8 *)(DataBuf + 11 ));
	volatile u8 aiTaskNumber =  (*(u8 *)(DataBuf + 12));
	volatile u8 aiExternalSyncEn = (*(u8 *)(DataBuf + 13));
	volatile u8 pulse_meas_en = 1;
	// if(aiExternalSyncEn&0x4) pulse_meas_en = 1;

	//volatile u8 aiTimeoutExSync = (*(u8 *)(DataBuf + 14));
	//volatile u8 RxGainFD_Array[3] = {255, 255, 255};
	//RxGainFD_Array[0] = (*(u8 *)(DataBuf + 14));
	//TimeoutDMA_mks = aiTimeoutExSync * 100000;


//	 xil_printf("splencode %d \n\r",splencode);
//	 xil_printf("aFDcount %d \n\r",aFDcount);
//	 xil_printf("aDetectThreshold %d \n\r",aDetectThreshold);
//	 xil_printf("aDelayCycle %d \n\r",aDelayCycle);
//	 xil_printf("aTimerPeriod %d \n\r",aTimerPeriod);
//	 xil_printf("aiAntCount %d \n\r",aiAntCount);
//	 xil_printf("aiTaskNumber %d \n\r",aiTaskNumber);
//	 xil_printf("aiExternalSyncEn %d \n\r",aiExternalSyncEn);
//	xil_printf("aTimeoutExSync %d \n\r",aTimeoutExSync);


		//volatile u16 *aiFmin = (u16*) malloc(aFDcount*sizeof(u16));
		//volatile u8 *aiFCount = (u8*) malloc(aFDcount*sizeof(u8));

		//if( (aiFmin ==NULL ) && (aiFCount == NULL) ){
			//xil_printf("Cmd %d malloc error!\n\r",CmdNum);
		//}

		volatile u16 i = 0;
		volatile u16 offset = 14;
	   for(i = 0;i<aFDcount;i++)
	   {

		   Pelengator.aiFmin[i] = (*(u16 *)(DataBuf + offset+0 + i*5)); // DataBuf[offset+1 + i*5] | (DataBuf[offset + i*4] << 8);
		   //Pelengator.aiFmin[i] = ReCalcFreq(Pelengator.aiFmin[i]);
		   Pelengator.aiFCount[i] = (*(u16 *)(DataBuf + offset+2 + i*5)); // DataBuf[offset+3 + i*5] | (DataBuf[offset+2 + i*4] << 8);
		   Pelengator.RxGain_Array[i] = (*(u8 *)(DataBuf + offset+4 + i*5)); //
		   //xil_printf("aiFmin[%d] = %d \n\r", i, Pelengator.aiFmin[i]);
		   //xil_printf("aiFCount[%d] = %d \n\r", i, Pelengator.aiFCount[i]);
		   //xil_printf("RxGain_Array[%d] = %d \n\r", i, Pelengator.RxGain_Array[i]);
	   }

	   volatile u8 TaskRepeatCount = (*(u8 *)(DataBuf + offset + aFDcount*(5)));
	   //xil_printf("TaskRepeatCount = %d \r\n", TaskRepeatCount);


	   volatile u8 FreqRepeatCount = (*(u8 *)(DataBuf + offset + aFDcount*(5) + 1));


	   //Адрес в пакете для количества повторений накопления
	   volatile u8 MaxHoldRepeatCount = (*(u8 *)(DataBuf + offset + aFDcount*(5) + 2));
	   volatile u8 locMaxHoldRepeatCount = MaxHoldRepeatCount;
	   //Адрес в пакете для задержки между циклами работы накопления в мс
	   volatile u16 MaxHoldDelay = (*(u16 *)(DataBuf + offset + aFDcount*(5) + 3));
	   //Введу заодно для аппаратного порога частоту и дельту
	   volatile u16 DetectFreq = (*(u16 *)(DataBuf + offset + aFDcount*(5) + 5));
	   volatile u16 DetectFreqDelta = (*(u16 *)(DataBuf + offset + aFDcount*(5) + 7));

	   Pelengator.DetectFreq = DetectFreq;
	   Pelengator.DetectFreqDelta = DetectFreqDelta;
	   //xil_printf("MaxHoldRepeatCount = %d \r\n", MaxHoldRepeatCount);
	   //xil_printf("MaxHoldDelay = %d \r\n", MaxHoldDelay);
	   //xil_printf("DetectFreq = %d \r\n", DetectFreq);
	   //xil_printf("DetectFreqDelta = %d \r\n", DetectFreqDelta);


	   volatile u8 locFreqRepeatCount = FreqRepeatCount;
	   if(FreqRepeatCount == 0)
	   {
		   FreqRepeatCount = 1;
		   //xil_printf("FreqRepeatCount def 0 = %d \r\n", FreqRepeatCount);
	   }
	   //else xil_printf("FreqRepeatCount = %d \r\n", FreqRepeatCount);
	   Pelengator.Len = splencode;



	   ////--------------
	   CmdNum = ETH_CMD_WB_SCAN;
	   //volatile u32 FreqStepNumber = Pelengator.FreqCount; //  116;
	   volatile u32 DataSize = 0; // (Pelengator.Len << 3) * FreqStepNumber ;
	   volatile u32 BytesSentOverall = 0;
	   volatile int idx;




	   volatile u8 FDcountINDEX = 0, ANTcountINDEX = 0, DetectEn = 0;
	   volatile u32 BytesSentReal = 0, DetectRezultFirstAnt = 0;
	   volatile s32 DetectRezult = 0;
	   volatile u16 AttenFreqIndex = 0;
	   volatile u16 AntAttenFreqIndex = 0;
	   volatile u8 local_rx_gain = 0;
	   //volatile u8 start_scan_delay = 0;
	   Pelengator.extra_byte_fft = 0;

	   u8 bfft_en = 1;
	   if(bfft_en)
	   {
		   WriteWord(rADC_DMA_MUX,2);//0 adc. 1 fft 2 fft complex
		   //Pelengator.extra_byte_fft = 4;
	   }
	   else
		{
		   WriteWord(rADC_DMA_MUX,0);
		   Pelengator.extra_byte_fft = 0;
		}



	   DataSize = 0;
	   u16 NStep = 0;

	   for(FDcountINDEX=0;FDcountINDEX < aFDcount; FDcountINDEX++)
	   {
		   DataSize += ((Pelengator.Len << Pelengator.ch_num_shift)+Pelengator.extra_byte_fft) * Pelengator.aiFCount[FDcountINDEX] * FreqRepeatCount  * MaxHoldRepeatCount * aiAntCount;
		   NStep = Pelengator.aiFCount[FDcountINDEX];
	   }
	   //volatile u32 locDataSize = DataSize;
	   //xil_printf("DataSize=%d \r\n", DataSize);
	   //Xil_DCacheFlushRange((UINTPTR)RX_BUFFER_BASE, DataSize);
	   Xil_DCacheInvalidateRange((UINTPTR)(RX_BUFFER_BASE),  DataSize);


	   volatile u8 rxgain = Pelengator.RxGain_Array[0]; // aiTimeoutExSync; // aiAntCount;




		if(rxgain != rxgain_prev)
		{
			rxgain_prev = rxgain;
			TALISE_setRxManualGain(&tal[TALISE_A], TAL_RX1,rxgain);
			TALISE_setRxManualGain(&tal[TALISE_A], TAL_RX2,rxgain);
			if(tal[TALISE_B].devHalInfo != NULL)
			{
				TALISE_setRxManualGain(&tal[TALISE_B], TAL_RX1,rxgain);
				TALISE_setRxManualGain(&tal[TALISE_B], TAL_RX2,rxgain);
			}
		}


	   if(lo_freq_hz_next != Pelengator.aiFmin[0]*1e6)
	   {
		   PrepareToGPIOHop(ReCalcFreq((u64)(Pelengator.aiFmin[0]*1e6)));
		   mdelay(1);
		   GPIO_TriggerHop();

		   mdelay(10);
		   GPIO_TriggerHop();
		   lo_freq_hz_next = Pelengator.aiFmin[0]*1e6;
	   }





	   while(TaskRepeatCount--)
	   {



		   BytesSentReal = 0;
		   BytesSentOverall = 0;
		   AttenFreqIndex = 0;
		   AntAttenFreqIndex = 0;
		   global_ind = AttenFreqIndex;



		   //volatile XTime t1, t2;
			//		   XTime_GetTime(&t1);

		   for(FDcountINDEX=0;FDcountINDEX < aFDcount; FDcountINDEX++)
		   {

			   //lo_freq_hz=  Pelengator.aiFmin[FDcountINDEX];
			   //lo_freq_hz *= 1000000;

			   rxgain = Pelengator.RxGain_Array[FDcountINDEX];
			   if(rxgain != rxgain_prev)
			   {
				   rxgain_prev = rxgain;
				   //TALISE_setRxManualGain(&tal[TALISE_A], TAL_RX1,rxgain);
				   //TALISE_setRxManualGain(&tal[TALISE_A], TAL_RX2,rxgain);
				   if(tal[TALISE_B].devHalInfo != NULL)
				   {
					   //TALISE_setRxManualGain(&tal[TALISE_B], TAL_RX1,rxgain);
				   	   //TALISE_setRxManualGain(&tal[TALISE_B], TAL_RX2,rxgain);
				   }
			   }




			   for(i=0;i < Pelengator.aiFCount[FDcountINDEX]; i++)
			   {
				   //LargeFreqStep(lo_freq_hz, rxgain);
				   //lo_freq_hz_last = lo_freq_hz;
				   lo_freq_hz = lo_freq_hz_next;

					   if( i == (Pelengator.aiFCount[FDcountINDEX]-1) )
				   {
					   if( FDcountINDEX == (aFDcount-1) )
					   {
						   lo_freq_hz_next = Pelengator.aiFmin[0]*1e6;
						  // xil_printf("FDcountINDEX == (aFDcount-1) \r\n");
					   }
					   else
					   {
						   lo_freq_hz_next = Pelengator.aiFmin[FDcountINDEX+1]*1e6;
						//   xil_printf("FDcountINDEX != (aFDcount-1) \r\n");
					   }
				   }
				   else
				   {
					   lo_freq_hz_next = lo_freq_hz+Pelengator.CalibrationFreqStep;
					//   xil_printf("i != (Pelengator.aiFCount[FDcountINDEX]-1) \r\n");
				   }

				   if(lo_freq_hz != lo_freq_hz_next)
				   {

					   PrepareToGPIOHop(ReCalcFreq(lo_freq_hz_next));//
					   GPIO_TriggerHop();

				   }

				   //printf("while lo_freq_hz = %d \r\n", (u32)((lo_freq_hz+Pelengator.CalibrationFreqStep)/1e6));

/*				uint64_t fhmRfPllFrequency_Hz;
					TALISE_getFhmRfPllFrequency(&tal[TALISE_A], &fhmRfPllFrequency_Hz);
					xil_printf("fhmRfPllFrequency_Hz = %d \r\n", (uint32_t)(fhmRfPllFrequency_Hz/1e6));
					xil_printf("lo_freq_hz = %d \r\n", (uint32_t)(lo_freq_hz/1e6));
					xil_printf("lo_freq_hz_next = %d \r\n", (uint32_t)(lo_freq_hz_next/1e6));
*/
				   locFreqRepeatCount = FreqRepeatCount;


				   if(pulse_meas_en) RstPulseMeas();

				   while(locFreqRepeatCount--)
				   {
					   global_FreqRepeatCount = locFreqRepeatCount;


					   //xil_printf("locFreqRepeatCount=%d\r\n", locFreqRepeatCount);
					   locMaxHoldRepeatCount = MaxHoldRepeatCount;
					   udelay(MaxHoldDelay);



					   while(locMaxHoldRepeatCount--)
					   {
						   global_MaxHoldRepeatCount = locMaxHoldRepeatCount;
						   //if(locMaxHoldRepeatCount%10 == 0) xil_printf("locMaxHoldRepeatCount=%d\r\n", locMaxHoldRepeatCount);


						   for(ANTcountINDEX=0;ANTcountINDEX < aiAntCount; ANTcountINDEX++)
						   {

							   global_AntCount=ANTcountINDEX;
							   SetExternalPorts(AttenFreqIndex, ANTcountINDEX);
							   udelay(aDelayCycle/**1e3*/);
							   //StartTimerDetect(((u32)aTimerPeriod)*1);
								// DelayMs(1);
							   //if( (locMaxHoldRepeatCount==0) && (locFreqRepeatCount==0) )
								//   xil_printf("00\r\n");
							   //while(1)
							   {


								   	if(bfft_en)
								   	{
										FFT_Rst();
										FFT_start(0);
								   	}


									DetectRezult = RecordAdcData(CmdNum, BytesSentOverall, bfft_en);
									Pelengator.ch0_blk_exp[AntAttenFreqIndex] = Pelengator.fftch4_status & 0xf;//AntAttenFreqIndex
									Pelengator.ch1_blk_exp[AntAttenFreqIndex] = (Pelengator.fftch4_status >> 8) & 0xf;//AntAttenFreqIndex
									//Pelengator.ch2_blk_exp[AttenFreqIndex] = (Pelengator.fftch4_status >> 16) & 0xf;
									//Pelengator.ch3_blk_exp[AttenFreqIndex] = (Pelengator.fftch4_status >> 24) & 0xf;
									//xil_printf("RecordAdcData\r\n");

									//RePack((UINTPTR)(RX_BUFFER_BASE+BytesSentOverall), AttenFreqIndex, AntAttenFreqIndex);//AntAttenFreqIndex
									//xil_printf("ProcessAdcDataTable 111\r\n");
									//FindSpectrumMax(RX_BUFFER_BASE+BytesSentOverall, ANTcountINDEX, AntAttenFreqIndex);

									//if(pulse_meas_en) GetPulseMeas(AttenFreqIndex);
									//FindSpectrumMax_PP(ANTcountINDEX);

									//ProcessAdcDataTable(RX_BUFFER_BASE+BytesSentOverall, Pelengator.Len, ANTcountINDEX);
									//FindSpectrumMax_test(RX_BUFFER_BASE+BytesSentOverall, ANTcountINDEX, AntAttenFreqIndex);
									//xil_printf("FindSpectrumMax AntAttenFreqIndex =%d ANTcountINDEX=%d\r\n", AntAttenFreqIndex, ANTcountINDEX);

							   }

							   BytesSentOverall = BytesSentOverall + (Pelengator.Len << Pelengator.ch_num_shift);// + Pelengator.extra_byte_fft;
							   if (DetectRezult == 0)
								   BytesSentReal = BytesSentReal + (Pelengator.Len << Pelengator.ch_num_shift) + Pelengator.extra_byte_fft;

							   AntAttenFreqIndex++;
						   }//for(ANTcountINDEX=0;ANTcountINDEX < aiAntCount; ANTcountINDEX++)


					   }// while(locMaxHoldRepeatCount--)



				   }//while(locFreqRepeatCount--)

				   if(pulse_meas_en) GetPulseMeas(AttenFreqIndex);
				   AttenFreqIndex++;
				   global_ind = AttenFreqIndex;

			   }//for(i=0;i < Pelengator.aiFCount[FDcountINDEX]; i++)



		   }//for(FDcountINDEX=0;FDcountINDEX < aFDcount; FDcountINDEX++)


		   //lo_freq_hz= Pelengator.CalibrationFreqInit;
		   //xil_printf("BytesSentReal = %d \r\n", BytesSentReal);


	   }//while(TaskRepeatCount--)



	   SetDirectExternalPorts();

	   //RegDiscreteClear(rLEDS,2);
		//XGpioPs_WritePin(&Gpio, 10, 0x0);
		//XGpioPs_WritePin(&Gpio, 12, 0x0);
		//DiscreteClr(rCPU_LEDs, 1);


	return 0;
}
void RePack(u8* buf, u16 freq_ind, u16 ant_freq_ind)
{
	//u64 t1 = TimemkS();
	memset(&fftpsd_rx1, 0, sizeof(fftpsd_rx1));
	memset(&fftpsd_rx2, 0, sizeof(fftpsd_rx2));

	u16 SpLen4 = Pelengator.Len/2;
	u16 SpLen2 = Pelengator.Len;
	u8 *address = (u8*)buf;
	s16 rx1_re, rx1_im, rx2_re, rx2_im;
	u16 j=0,i=0;
	double f_rx1_re, f_rx1_im, f_rx1_mag, f_rx1_mag_log, f_rx1_mag_max=0;
	double f_rx2_re, f_rx2_im, f_rx2_mag, f_rx2_mag_log, f_rx2_mag_max=0;
	s64 s64_rx1_re, s64_rx1_im, s64_rx2_re, s64_rx2_im;
	comp tmpcomp_rx1;
	comp tmpcomp_rx2;

    u8 ind_ch0_I = 0; // 0;
    u8 ind_ch0_Q = 2; // 2;
    u8 ind_ch1_I = 4; // 4;
    u8 ind_ch1_Q = 6; // 6;



    for (j=0; j<SpLen4; j++)
    {
		rx1_re = *(s16*) (address + j*8 + ind_ch0_I);
		rx1_im = *(s16*) (address + j*8 + ind_ch0_Q);




		s64_rx1_re = (s64)rx1_re;
		s64_rx1_im = (s64)rx1_im;
		s64_rx1_re = s64_rx1_re << Pelengator.ch0_blk_exp[ant_freq_ind];
		s64_rx1_im = s64_rx1_im << Pelengator.ch0_blk_exp[ant_freq_ind];
		f_rx1_re = (double) s64_rx1_re;
		f_rx1_im = (double) s64_rx1_im;
		//f_rx1_re *= (pow(2,Pelengator.ch0_blk_exp[ant_freq_ind]));
		//f_rx1_im *= (pow(2,Pelengator.ch0_blk_exp[ant_freq_ind]));

//		f_rx1_mag = sqrt(f_rx1_re*f_rx1_re + f_rx1_im*f_rx1_im);
//		if(f_rx1_mag != 0)
//			f_rx1_mag_log = 20*log10(f_rx1_mag)-100-30;
//		else
//			f_rx1_mag_log = -240;


		rx2_re = *(s16*) (address + j*8 + ind_ch1_I);
		rx2_im = *(s16*) (address + j*8 + ind_ch1_Q);

		s64_rx2_re = (s64)rx2_re;
		s64_rx2_im = (s64)rx2_im;
		s64_rx2_re = s64_rx2_re << Pelengator.ch1_blk_exp[ant_freq_ind];
		s64_rx2_im = s64_rx2_im << Pelengator.ch1_blk_exp[ant_freq_ind];


		f_rx2_re = (double)	s64_rx2_re;
		f_rx2_im = (double) s64_rx2_im;

		//f_rx2_re *= (pow(2,Pelengator.ch1_blk_exp[ant_freq_ind]));
		//f_rx2_im *= (pow(2,Pelengator.ch1_blk_exp[ant_freq_ind]));


//		f_rx2_mag = sqrt(f_rx2_re*f_rx2_re + f_rx2_im*f_rx2_im);
//		if(f_rx2_mag != 0)
//			f_rx2_mag_log = 20*log10(f_rx2_mag)-100-30;
//		else
//			f_rx2_mag_log = -240;


		tmpcomp_rx1 = f_rx1_re + I*f_rx1_im;
		tmpcomp_rx2 = f_rx2_re + I*f_rx2_im;
       if (Pelengator.RXcAA[freq_ind] != 0)
       {
    	   fftpsd_rx1[SpLen4-j] = conj(tmpcomp_rx1);
    	   fftpsd_rx2[SpLen4-j] = conj(tmpcomp_rx2);
       }
       else
       {
    	   fftpsd_rx1[j + SpLen4] = tmpcomp_rx1;
    	   fftpsd_rx2[j + SpLen4] = tmpcomp_rx2;
       }
       i++;
       //xil_printf("f_rx1_re=%d,  f_rx1_im=%d   ",f_rx1_re,f_rx1_im);////////////////////////////////////////////////////////////////////////////////////////////

    };

    for (j=SpLen4; j<SpLen2; j++)
    {

		rx1_re = *(s16*) (address + j*8 + ind_ch0_I);
		rx1_im = *(s16*) (address + j*8 + ind_ch0_Q);

		s64_rx1_re = (s64)rx1_re;
		s64_rx1_im = (s64)rx1_im;
		s64_rx1_re = s64_rx1_re << Pelengator.ch0_blk_exp[ant_freq_ind];
		s64_rx1_im = s64_rx1_im << Pelengator.ch0_blk_exp[ant_freq_ind];

		f_rx1_re = (double) s64_rx1_re;
		f_rx1_im = (double) s64_rx1_im;
		//f_rx1_re *= (pow(2,Pelengator.ch0_blk_exp[ant_freq_ind]));
		//f_rx1_im *= (pow(2,Pelengator.ch0_blk_exp[ant_freq_ind]));

		//f_rx1_mag = sqrt(f_rx1_re*f_rx1_re + f_rx1_im*f_rx1_im);
//		if(f_rx1_mag != 0)
//			f_rx1_mag_log = 20*log10(f_rx1_mag)-100-30;
//		else
//			f_rx1_mag_log = -240;


		rx2_re = *(s16*) (address + j*8 + ind_ch1_I);
		rx2_im = *(s16*) (address + j*8 + ind_ch1_Q);


		s64_rx2_re = (s64)rx2_re;
		s64_rx2_im = (s64)rx2_im;
		s64_rx2_re = s64_rx2_re << Pelengator.ch1_blk_exp[ant_freq_ind];
		s64_rx2_im = s64_rx2_im << Pelengator.ch1_blk_exp[ant_freq_ind];

		f_rx2_re = (double) s64_rx2_re;
		f_rx2_im = (double) s64_rx2_im;

		//f_rx2_re *= (pow(2,Pelengator.ch1_blk_exp[ant_freq_ind]));
		//f_rx2_im *= (pow(2,Pelengator.ch1_blk_exp[ant_freq_ind]));


//		f_rx2_mag = sqrt(f_rx2_re*f_rx2_re + f_rx2_im*f_rx2_im);
//		if(f_rx2_mag != 0)
//			f_rx2_mag_log = 20*log10(f_rx2_mag)-100-30;
//		else
//			f_rx2_mag_log = -240;


		tmpcomp_rx1 = f_rx1_re + I*f_rx1_im;
		tmpcomp_rx2 = f_rx2_re + I*f_rx2_im;




       if (Pelengator.RXcAA[freq_ind] != 0)
       {
     	   fftpsd_rx1[SpLen2 + SpLen4 - j] = conj(tmpcomp_rx1);
     	   fftpsd_rx2[SpLen2 + SpLen4 - j] = conj(tmpcomp_rx2);

       }
       else
       {
     	   fftpsd_rx1[j - SpLen4] = tmpcomp_rx1;
     	   fftpsd_rx2[j - SpLen4] = tmpcomp_rx2;
       }

    };

  //  xil_printf("time repack = %d \r\n", (u32)((TimemkS() - t1)));

}

void apply_blackman_window(s16 *input_win, int length) {

    for (int n = 0; n < length; n++) {
    	//input2[n]=input_win[n];

    	input_win[n] *= 0.42 - 0.5 * cos(2 * M_PI * n / (length - 1)) + 0.08 * cos(4 * M_PI * n / (length - 1));
    }
}







void compute_phase_unwrap(s16 *I_fase, s16 *Q_fase, double *phase_unwrapped_deg, u32 *jump_indices, u32 *num_jumps, u32 length, u32 threshold) {
    const double RAD_TO_DEG = 180.0 / M_PI;
    double prev_phase_rad = 0;
    double phase_rad, phase_diff_rad;

    *num_jumps = 0;
    u16 koef_accum=1;
    for (int i = 0; i < length; i++) {
        // Вычисление фазы в радианах
        phase_rad = atan2(Q_fase[i], I_fase[i]);

        char buffer[50];
        sprintf(buffer, "%f", phase_rad);
		xil_printf("phase_rad[%d] = %s \n", i, buffer );

		usleep(5000);




        // Применение фазовой развёртки
        if (i > 0) {
            phase_diff_rad = phase_rad - prev_phase_rad;
            if (phase_diff_rad > M_PI) {
                phase_rad -= 2 * M_PI;
            } else if (phase_diff_rad < -M_PI) {
                phase_rad += 2 * M_PI*koef_accum;
                if (phase_rad<prev_phase_rad)
                {
                	koef_accum++;
                	phase_rad += 2 * M_PI;
                }
            }
        }

        // Определение и сохранение точек скачков фазы
        if (i > 0) {
            phase_diff_rad = phase_rad - prev_phase_rad;
            // Игнорируем скачки, вызванные циклическим переходом от 360° к 0°
            if (fabs(phase_diff_rad) > threshold && fabs(phase_diff_rad) < 2 * M_PI) {
                jump_indices[*num_jumps] = i;
                (*num_jumps)++;
            }
        }
        // Перевод развёрнутой фазы в градусы
        phase_unwrapped_deg[i] = phase_rad * RAD_TO_DEG;
        prev_phase_rad = phase_rad;
    }
}

#define N2 4  // Порядок фильтра (можно изменить)
#define SCALE_FACTOR 32768  // Для работы с целыми числами (если нужен фиксированный формат)

// Коэффициенты фильтра (пример для НЧ-фильтра Баттерворта 4-го порядка, нормализованного)
int16_t a[4+1] = {32768, -31855, 18290, -6040, 722};  // Обратные коэффициенты
int16_t b[4+1] = {585, 2342, 3514, 2342, 585};  // Прямые коэффициенты

// Буферы задержек
int32_t x[N2+1] = {0};  // Входные отсчёты
int32_t y[N2+1] = {0};  // Выходные отсчёты

// Функция обработки сигнала через фильтр
s16 butterworth_filter(s16 input) {
    int32_t acc = 0;

    // Сдвиг буферов
    for (int i = N2; i > 0; i--) {
        x[i] = x[i-1];
        y[i] = y[i-1];
    }
    x[0] = input;

    // Применение разностного уравнения
    for (int i = 0; i <= N2; i++) {
        acc += (b[i] * x[i]) - (a[i] * y[i]);
    }

    // Масштабирование результата
    y[0] = acc / SCALE_FACTOR;

	char str1[20];
	sprintf(str1, "%.6d", input);
	printf("signal = %s,\n",str1);
    //
	 str1[20]=0;
	sprintf(str1, "%.12d", y[0]);
	printf("output = %s,\n",str1);
    return (int16_t)y[0];
}


void read_ddr3_reports(u32 base_addr) {//RX_BUFFER_BASE

	u32 count_o = N1;//65536
	Xil_DCacheInvalidateRange((UINTPTR)(base_addr),  count_o*16);
	memset(&BUF_ADC_I, 0, sizeof(BUF_ADC_I));//adcpsd_rx1
	memset(&BUF_ADC_Q, 0, sizeof(BUF_ADC_Q));//adcpsd_rx1
	memset(&BUF_ADC_I_filter, 0, sizeof(BUF_ADC_I_filter));//adcpsd_rx1
	memset(&BUF_ADC_Q_filter, 0, sizeof(BUF_ADC_Q_filter));//adcpsd_rx1

	u32 index1=base_addr;
	s16 data2_Q;
	s16 data2_I;

	u32 k=0;

    while(k<count_o) {
    	data2_Q = (s16) Xil_In16(index1+16*k+6);
    	data2_I = (s16) Xil_In16(index1+16*k+4);
//    	if (k<20){
//    		xil_printf("dataI %d \n",data2_I);
//    	}
    	BUF_ADC_Q[k]=  ((float )(data2_Q));   //*2.5/5345);//Какая частота дескритизации выборок?245
    	BUF_ADC_I[k]= ((float )(data2_I));     //*2.5/5345);


    	k+=1;
    }


    xil_printf("CHEK MODULATIONS \n");
//	CHEK_FASE_MOD(BUF_ADC_I, BUF_ADC_Q, 10e6, (245000000/4));
//	CHEK_LORA_MOD(BUF_ADC_I, BUF_ADC_Q, 3e6, (245000000/4));
	CHEK_LORA_MOD_2(BUF_ADC_I, BUF_ADC_Q, 3e6, (245000000/4));

}

//////////конец фукнкции!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//int save_buffer_to_file(const char *filename, uint32_t *buffer, uint32_t size) {
//	xil_printf("file is create \n");
//	FILE *file;
//
//	// Открытие файла для записи (режим "w" создаёт файл, если он не существует)
//	file = fopen("filename", "w");
//
//
//
//    fclose(file);
//    xil_printf("close file \n");
//}







void FindSpectrumMax_test(u8 *buf, u8 AntNum, u8 ind)
{
	double f_rx1_re, f_rx1_im, f_rx1_mag, f_rx1_mag_log;
	double f_rx2_re, f_rx2_im, f_rx2_mag, f_rx2_mag_log;

	double sample=-240, sample_rx2=-240, sample_max=-240;
	u16 ind_max=0;

	 comp rx1_comp ;
	 comp rx2_comp ;
	volatile u16 i=0;
	for (i = Pelengator.StartAnalyse; i < Pelengator.StopAnalyse; i++)
	//for (i = 0; i < 256; i++)
    {
		f_rx1_mag = cabs(fftpsd_rx1[i]); // sqrt(f_rx1_re*f_rx1_re + f_rx1_im*f_rx1_im);
		if(f_rx1_mag != 0)
			f_rx1_mag_log = 20*log10(f_rx1_mag)-100-30;
		else
			f_rx1_mag_log = -240;

		f_rx2_mag = cabs(fftpsd_rx2[i]); // sqrt(f_rx2_re*f_rx2_re + f_rx2_im*f_rx2_im);
		if(f_rx2_mag != 0)
			f_rx2_mag_log = 20*log10(f_rx2_mag)-100-30;
		else
			f_rx2_mag_log = -240;


		sample = f_rx1_mag_log;
		sample_rx2 = f_rx2_mag_log;
		if(sample > sample_max)
		{
			sample_max = sample;
			ind_max = i;
			rx1_comp = f_rx1_re + I*f_rx1_im;
			rx2_comp = f_rx2_re + I*f_rx2_im;
		}
		if(i==246)
			xil_printf("246 freq = %d power1 = %d power2 = %d \r\n", (s32)(FFTcodeToFreq_1(i)+lo_freq_hz/1e6), (s32)sample, (s32)sample_rx2);
		//xil_printf("i=%d freq = %d power1 = %d power2 = %d \r\n",i, (s32)(FFTcodeToFreq_1(i)+lo_freq_hz/1e6), (s32)sample, (s32)sample_rx2);
    }



	double freqtmp = FFTcodeToFreq_1(ind_max)+lo_freq_hz/1e6;
	//printf("AddSignal freq = %.2f power1 = %.2f power2 = %.2f antnum=%d \r", freqtmp, sample, sample_rx2, AntNum);
	xil_printf("AddSignal freq = %d power1 = %d power2 = %d antnum=%d i=%d \r\n", (s32)freqtmp, (s32)sample_max, (s32)sample_rx2, AntNum, i);
	//xil_printf("re = %d im = %d re2 = %d im2=%d ", rx1_re, rx1_im, rx2_re, rx2_im);
	//xil_printf("ch0_blk_exp = %d ch1_blk_exp = %d  \r\n", Pelengator.ch0_blk_exp[ind], Pelengator.ch1_blk_exp[ind]);
    if(sample_max > Pelengator.MaxThreshold)
	AddSignal(
        freqtmp,
        sample,
		sample_rx2,
		rx2_comp,
		rx1_comp,
        AntNum
        );


}
void FindSpectrumMax(u8 *buf, u8 AntNum, u8 ind)
{
	//xil_printf("FindSpectrumMax antnum=%d \r\n", AntNum);
	//if(lo_freq_hz==2490e6) if(AntNum==0)printf("FindSpectrumMax \r");
	//u8 *address = (u8*)buf;
	//s16 rx1_re, rx1_im, rx2_re, rx2_im;
	double f_rx1_re, f_rx1_im, f_rx1_mag, f_rx1_mag_log;
	double f_rx2_re, f_rx2_im, f_rx2_mag, f_rx2_mag_log;
	//comp rx1_comp;
	//comp rx2_comp;
	double sample_rx2=0;

	float GlobalThreshold = Pelengator.MaxThreshold; // -55.0;//StrToFloat(formParameters->leThreshold->Text);
    int TentsCeiling = Pelengator.TensCeiling; // 20;//StrToInt(formParameters->leTentsCeiling->Text);
    int AbsenceLimit = Pelengator.Absence; // 30;//StrToInt(formParameters->leAbsence->Text);
    int ShiftC = 0;//GetPrmInt(3);
    int StartAnalyse = Pelengator.StartAnalyse ;// 0; // 0;
    int StopAnalyse = Pelengator.StopAnalyse; // round(Pelengator.Len*FreqBand/Fd); // SamplesNumber;///UT.iAntCount;

	u16 SamplesNumber = Pelengator.Len;
	u8 SignalTentBuffer[SamplesNumber];
    memset(SignalTentBuffer,0,SamplesNumber);
    double sample=0;
    double ampj = 0;
    int indj = 0;
    int ci = 1;//int
    double cilog=0;
    int fi;
    int AbsenceCounter=0;

   // printf("GlobalThreshold =  %.2f \r", GlobalThreshold);
	volatile u16 i=0, k = 0;
	for (i = StartAnalyse; i < StopAnalyse; i++)
    {
		//rx1_re = *(s16*) (address + i*8 + 0);
		//rx1_im = *(s16*) (address + i*8 + 2);
        //
		//f_rx1_re = (double) rx1_re;
		//f_rx1_im = (double) rx1_im;
		//f_rx1_re *= (pow(2,Pelengator.ch0_blk_exp[ind]));
		//f_rx1_im *= (pow(2,Pelengator.ch0_blk_exp[ind]));

		f_rx1_mag = cabs(fftpsd_rx1[i]); // sqrt(f_rx1_re*f_rx1_re + f_rx1_im*f_rx1_im);
		if(f_rx1_mag != 0)
			f_rx1_mag_log = 20*log10(f_rx1_mag)-100-30;
		else
			f_rx1_mag_log = -240;


		//rx2_re = *(s16*) (address + i*8 + 4);
		//rx2_im = *(s16*) (address + i*8 + 6);
        //
		//f_rx2_re = (double) rx2_re;
		//f_rx2_im = (double) rx2_im;
        //
		//f_rx2_re *= (pow(2,Pelengator.ch1_blk_exp[ind]));
		//f_rx2_im *= (pow(2,Pelengator.ch1_blk_exp[ind]));


		f_rx2_mag = cabs(fftpsd_rx2[i]); // sqrt(f_rx2_re*f_rx2_re + f_rx2_im*f_rx2_im);
		if(f_rx2_mag != 0)
			f_rx2_mag_log = 20*log10(f_rx2_mag)-100-30;
		else
			f_rx2_mag_log = -240;


		sample = f_rx2_mag_log; // SamplesBuf[i];

		//if(i==10)printf("sample[%d] =  %.2f \r",i, sample);
		fi = i-indj;
		if (fi<TentsCeiling) ci=1;
		else ci=fi/TentsCeiling;

		cilog= 20*log10(ci);

		if (sample<=GlobalThreshold)
		{
			if (AbsenceCounter++==AbsenceLimit) ampj=-100;
		}
		else AbsenceCounter=0;
		//xil_printf("AbsenceCounter[%d] = %d\r\n",i, AbsenceCounter);
		//if (sample*ci>=ampj) {
		if (sample+cilog>=ampj)
		{
			indj=i-ShiftC;
			ampj=sample;
			SignalTentBuffer[i] = 1;
		}
		else SignalTentBuffer[i] = 0;
    } // for (i = StartAnalyse; i < StopAnalyse; i++)

	ampj = 0;
    indj = SamplesNumber-1;
    AbsenceCounter = 0;
    ci = 1;
   // printf("(i = StopAnalyse; i >= StartAnalyse; i--)\r");
    //if(lo_freq_hz==2490e6) if(AntNum==0)xil_printf("StopAnalyse =  %d StartAnalyse= %d\r\n", StopAnalyse, StartAnalyse);
	for (i = StopAnalyse; i >= StartAnalyse; i--)
    {
		//rx1_re = *(s16*) (address + i*8 + 0);
		//rx1_im = *(s16*) (address + i*8 + 2);
        //
		//f_rx1_re = (double) rx1_re;
		//f_rx1_im = (double) rx1_im;
		//f_rx1_re *= (pow(2,Pelengator.ch0_blk_exp[ind]));
		//f_rx1_im *= (pow(2,Pelengator.ch0_blk_exp[ind]));

		f_rx1_mag = cabs(fftpsd_rx1[i]); //sqrt(f_rx1_re*f_rx1_re + f_rx1_im*f_rx1_im);
		if(f_rx1_mag != 0)
			f_rx1_mag_log = 20*log10(f_rx1_mag)-100-30;
		else
			f_rx1_mag_log = -240;

		//rx2_re = *(s16*) (address + i*8 + 4);
		//rx2_im = *(s16*) (address + i*8 + 6);
        //
		//f_rx2_re = (double) rx2_re;
		//f_rx2_im = (double) rx2_im;
        //
		//f_rx2_re *= (pow(2,Pelengator.ch1_blk_exp[ind]));
		//f_rx2_im *= (pow(2,Pelengator.ch1_blk_exp[ind]));

		f_rx2_mag = cabs(fftpsd_rx2[i]); //sqrt(f_rx2_re*f_rx2_re + f_rx2_im*f_rx2_im);
		if(f_rx2_mag != 0)
			f_rx2_mag_log = 20*log10(f_rx2_mag)-100-30;
		else
			f_rx2_mag_log = -240;




		sample = f_rx2_mag_log; // SamplesBuf[i];
		sample_rx2 = f_rx1_mag_log; // SamplesBuf[i];
		comp rx1_comp = f_rx1_re + I*f_rx1_im;
		comp rx2_comp = f_rx2_re + I*f_rx2_im;
		//if(lo_freq_hz==2490e6) if(AntNum==0) xil_printf("sample[%d] =  %d \r",i, (s32)sample);

		fi = indj-i;
		if (fi<TentsCeiling) ci=1;
		else ci=fi/TentsCeiling;

		cilog= 20*log10(ci);

		if (sample<=GlobalThreshold)
		{
			if (AbsenceCounter++==AbsenceLimit) ampj=-100;
		}
		else AbsenceCounter=0;
		//if (sample*ci>=ampj) {
		if (sample+cilog>=ampj)
		{
			if ((i>=StartAnalyse)&&(i<StopAnalyse)&&(sample>GlobalThreshold)&&(SignalTentBuffer[i]!=0))
			{

				double freqtmp = FFTcodeToFreq_1(i)+lo_freq_hz/1e6;
				//printf("AddSignal freq = %.2f power1 = %.2f power2 = %.2f antnum=%d \r", freqtmp, sample, sample_rx2, AntNum);
				//xil_printf("AddSignal freq = %d power1 = %d power2 = %d antnum=%d i=%d \r\n", (s32)freqtmp, (s32)sample, (s32)sample_rx2, AntNum, i);
				//xil_printf("re = %d im = %d re2 = %d im2=%d ", rx1_re, rx1_im, rx2_re, rx2_im);
				//xil_printf("ch0_blk_exp = %d ch1_blk_exp = %d  \r\n", Pelengator.ch0_blk_exp[ind], Pelengator.ch1_blk_exp[ind]);
                if(AddSignal(
                    freqtmp,
                    sample,
					sample_rx2,
					fftpsd_rx1[i], //rx2_comp,
					fftpsd_rx2[i], // rx1_comp,
                    AntNum
                    ) > 0)
                {
                    //Sp->SpectrumMarkerArray.Add(TMySpectrumMarker(FreqVector[i]));
					// добавлен сигнал в таблицу


                }



			}//if ((i>=StartAnalyse)&&(i<StopAnalyse)&&(sample>GlobalThreshold)&&(SignalTentBuffer[i]!=0))
			indj=i+ShiftC;
			ampj=sample;
		}
	}//for (i = StopAnalyse; i >= StartAnalyse; i--)

	return;
}




int AddSignal(double Freq, float Power, float Power_rx2, comp CompCh1, comp CompCh2, int AntNum)
{
	xil_printf("AddSignal AntNum=%d Freq=%d Power=%d\r\n", AntNum, (s32)Freq,(s32)Power);
	//if (CheckSpurs(Freq,Power)>0)
	{
		struct TSignal Signal;
		memset(&Signal, 0, sizeof(TSignal));
      //Signal = new sSignal;
      Signal.Freq = Freq;
     // if(AntNum==0)
      {
    	  Signal.Power = Power;
    	  Signal.power = (s8)round(Power);
      }
      Signal.Sch1[AntNum] = CompCh1;
      Signal.Sch2[AntNum] = CompCh2;
      Signal.AntNum[AntNum] = AntNum;
      Signal.PelengFi = 0;
      Signal.PelengTeta = 0;

      Signal.fc_min = (u16)round(Freq);
      Signal.az = 0;//-180
      Signal.el = 0;


      //xil_printf("Signal.power = %d \r\n", Signal.power);
      volatile u16 i=0;

      /*if(fabs(Freq - 4000.0) < 2)
      {
    	  xil_printf("Signal.power = %d ", (s32)Signal.power);
    	  double tmp = cabs(CompCh1);
    	  if(tmp != 0)
    		  tmp  = 20*log10(tmp)-100-30;
    	  	else
    	  		tmp = 0;

    	  xil_printf("addtmp[%d]=%d ", AntNum, (s32)tmp);
    	  tmp = cabs(CompCh2);
    	   if(tmp != 0)
    		   tmp  = 20*log10(tmp)-100-30;
    	   else
    		   tmp = 0;
    	   xil_printf(" %d \r\n",(s32)tmp);
      	}*/
	/*	double tmp = cabs(CompCh1);
		if(tmp != 0)
			tmp  = 20*log10(tmp)-100-30;
		else
			tmp = 0;

		xil_printf("tmp[%d]=%d \n\r", AntNum, (s32)tmp);

*/
      //Проверка наличия сигнала на данной частоте в списке сигналов
      for (i=0; i< Pelengator.SignalCnt;i++)
      {
    	  TSignal * TmpSignal=&(SignalTable[i]);//(sSignal*)SignalList->Items[i];
    	  float Tmp = fabs(Freq - TmpSignal->Freq);

    	  if (Tmp < 3)
    	  {
         	//Если сигнал на данной частоте в списке присутствует, то обновляем мощность и время


         	TmpSignal->DetectTime = TimeS();
         	//if(AntNum==0)
         	{
         		TmpSignal->Power = Power;
            	TmpSignal->power = (s8)round(Power);
         	}
            TmpSignal->imp_len_max = TmpSignal->DetectionCounter;
            TmpSignal->DetectionCounter++;
            TmpSignal->Sch1[AntNum] = Signal.Sch1[AntNum];
            TmpSignal->Sch2[AntNum] = Signal.Sch2[AntNum];
            TmpSignal->AntNum[AntNum] = AntNum;
            TmpSignal->Refresh=1;

            /*TmpSignal->imp_period = Signal->imp_period;
            TmpSignal->imp_duration = Signal->imp_duration;
            TmpSignal->imp_ampmax = Signal->imp_ampmax;*/

            //if (AntNum == 0) TmpSignal->AntNum = 0;
            //else TmpSignal->AntNum = TmpSignal->AntNum + 1;
            //RefreshSignalTable();
            //delete Signal;
            return 1;
    	  }
      }
      //if(fabs(Freq - 4000.0)<1)
      	//xil_printf("addnew 4=%d\r\n",(s32)Freq);
      //Сигнала на данной частоте нет => добавляем в список
      Signal.DetectTime = TimeS();
      Signal.DetectionCounter=1;
      Signal.imp_len_max = 1;
      //if(SignalList->Count<10)
          //SignalList->Add(Signal);
      //xil_printf("new AddSignal freq=%d power=%d\r\n",(u32)Freq, (s32)Power);
          AddToSignalTable(Signal);

      //RefreshSignalTable();
      //SignalThreadList->UnlockList();
      return 1;
   }
   //else {return 0;}
}
void AddToSignalTable(TSignal S)
{
	volatile u16 i=0;
	for(i=Pelengator.SignalCnt;i>0;i--)
	{
		SignalTable[i] = SignalTable[i-1];
		//SignalTable[i].imp_len_min = i;
	}
	SignalTable[0] = S;
	if(Pelengator.SignalCnt < 1000/*sizeof(SignalTable)*/) Pelengator.SignalCnt++;
}
void DeleteFromSignalTable(u16 ind)
{
	volatile u16 i=0;
	bool find=0;
	if (ind==0) find = 1;
	for(i=0;i<Pelengator.SignalCnt;i++)//-1//sizeof(SignalTable)
	{
		if(i==ind) find=1;
		if(find==1)
		{
			//xil_printf("del s[%d].freq =%d\r\n",i, (s32)SignalTable[i].Freq);
			SignalTable[i] = SignalTable[i+1];
		}
	}

	if(Pelengator.SignalCnt > 0) Pelengator.SignalCnt--;
	return;
}
void CheckSignalsTime()
{
        TSignal *Signal;
        //TList *SignalList = SignalThreadList->LockList();
        //double DeleteTime = 10.0/24/3600;
        double DeleteTime = Pelengator.AbsenceTime;///24/3600;     // 10.0/24/3600;//
        double CurrentTime = (double)TimeS()-DeleteTime;
        int StartIdx = Pelengator.SignalCnt-1;

        for (int i=StartIdx; i>-1; i--)
        //for (int i=0; i<count; i++)
        {
                Signal = &SignalTable[i];

                double TimeDif = CurrentTime - (double)(Signal->DetectTime);
                //xil_printf("Signal->DetectTime =%d\r\n", (s32)Signal->DetectTime);
                //xil_printf("CurrentTime =%d\r\n", (s32)CurrentTime);
                //xil_printf("TimeDif[%d] =%d\r\n",i, (s32)TimeDif);
                if (TimeDif > 0)
                {
                	DeleteFromSignalTable(i);

                }
        }
        //SignalThreadList->UnlockList();
}
u8 CheckSignalAnt(TSignal * S)
{
	//xil_printf("(*S).freq =%d\r\n", (s32)(*S).Freq);
	u8 err=0;
	for(int i=0;i<11;i++)
	{
		//xil_printf("(*S).AntNum[%d] =%d\r\n",i, (s32)(*S).AntNum[i]);
		if((*S).AntNum[i] != i)
			//return 0;
			err++;

	}
	if(err>0)// минимальное количество антенн 0 сигнал должен быть во всех антеннах, 1 - может на одной не быть
		return 0;
	else
		return 1;

}
u8 ClearSignalPeleng(TSignal * S)
{
	for(int i=0;i<11;i++)
	{
		(*S).AntNum[i]=0;
	}
	(*S).az=0;
	(*S).el=0;
	return 1;
}
void ClearSignalTable()
{
//xil_printf("ClearSignalTable cnt=%d\r\n", Pelengator.SignalCnt);
        for (int i=0; i<Pelengator.SignalCnt+1; i++)
        //for (int i=0; i<count; i++)
        {
        	//xil_printf("s[%d].freq =%d\r\n",i, (s32)SignalTable[i].Freq);
        	DeleteFromSignalTable(i);

        }
        //xil_printf("ClearSignalTable end cnt=%d\r\n", Pelengator.SignalCnt);
}
void FilterSignalTable()
{
//xil_printf("ClearSignalTable cnt=%d\r\n", Pelengator.SignalCnt);
        for (int i=0; i<Pelengator.SignalCnt+1; i++)
        //for (int i=0; i<count; i++)
        {
        	//xil_printf("s[%d].freq =%d\r\n",i, (s32)SignalTable[i].Freq);
        	//if(SignalTable[i].imp_len_min <= 1)
        	if(SignalTable[i].Power <= -40.0)
        		DeleteFromSignalTable(i);

        }
        //xil_printf("ClearSignalTable end cnt=%d\r\n", Pelengator.SignalCnt);
}
void CalculatePeleng(void)
{
	//u32 t1 = TimemS();

	//xil_printf("CalculatePeleng Pelengator.SignalCnt=%d\n\r", Pelengator.SignalCnt);
	double Um = 0.0;
	double peleng = 0.0;
	u8 antnum = 11;
	double tmp;
	comp S[11];
	comp O[11];

    int StartIdx = Pelengator.SignalCnt-1;

    for (int i=StartIdx; i>-1; i--)
	//for(int i=0;i<Pelengator.SignalCnt;i++)//sizeof(SignalTable)
	{
		//SignalTable[i] = SignalTable[i+1];

		//if(fabs(SignalTable[i].Freq - 4000.0) < 1)
		//if(i<10)
		{
			//xil_printf("SignalTable[%d] freq=%d power=%d AntNumCnt=%d \r\n", i, (u32)SignalTable[i].Freq, (s32)SignalTable[i].Power, SignalTable[i].AntNumCnt);


			peleng = -1;
			Um = -1;
			if(CheckSignalAnt(&SignalTable[i]))
			{
				for(int k=0;k<11;k++)
				{

					tmp = cabs(SignalTable[i].Sch1[k]);
					if(tmp != 0)
						tmp  = 20*log10(tmp)-100-30;
					else
						tmp = 0;

					//xil_printf("SignalTable[%d].abs Sch1[%d]=%d \n\r", i, k, (s32)tmp);

					tmp = cabs(SignalTable[i].Sch2[k]);
					if(tmp != 0)
						tmp  = 20*log10(tmp)-100-30;
					else
						tmp = 0;

					//xil_printf("SignalTable[%d].abs Sch2[%d]=%d \n\r", i, k, (s32)tmp);
					S[k] = SignalTable[i].Sch1[k];// основной от коммутатора
					O[k] = SignalTable[i].Sch2[k];// опорный ценральный
				}

				//xil_printf("i =%d CheckSignalAnt true \n\r",i);
				peleng = calc_peleng_car(SignalTable[i].Freq, S, O, &Um, antnum);
				SignalTable[i].az = (u16)(peleng*10+1800);//-180
				SignalTable[i].el = (u16)Um;
				//xil_printf("ind = %d freq = %d ; int az =%d el=%d \r\n", i, (u32)SignalTable[i].Freq, SignalTable[i].az, SignalTable[i].el);
				//xil_printf("double az =%d el=%d \r\n", (s32)peleng, (s32)Um);
				//ClearSignalPeleng(&SignalTable[i]);
			}
			else
			{
				//xil_printf("i =%d CheckSignalAnt false \n\r",i);
				DeleteFromSignalTable(i);
			}
			//xil_printf("\r\n");
			//xil_printf("az =%d el=%d \r\n ", (s32)peleng, (s32)Um);
			//xil_printf("az =%d el=%d \n\r", (s32)peleng, (s32)Um);
			//xil_printf("peleng = %d; Um =  %d \r\n",(s32)Pelengator.peleng_angle2, (s32)Pelengator.peleng_angle);
		}

	}
	// xil_printf("filter Pelengator.SignalCnt = %d \r\n", Pelengator.SignalCnt);
	/*comp vecS[11];
	comp vecO[11];
	double Um = 0.0;
	 double peleng = calc_peleng_car(FFTcodeToFreq(rx1_mag_max_ind) + Pelengator.aiFmin[0], vecS, vecO, &Um, antnum);
	 Pelengator.peleng_angle = Um;
	 Pelengator.peleng_angle2 = peleng;*/
   // xil_printf("time calc peleng = %d \r\n", (u32)(TimemS() - t1));
	return;
}
print_SignalTable()
{
	xil_printf("print_SignalTable Pelengator.SignalCnt = %d \r\n", Pelengator.SignalCnt);
	for(int i=0;i<10/*Pelengator.SignalCnt*/;i++)//sizeof(SignalTable)
	{
		//SignalTable[i] = SignalTable[i+1];
//		xil_printf("SignalTable[%d] freq=%d power=%d cnt=%d \n\r", i, (u32)SignalTable[i].Freq, (s32)SignalTable[i].power, SignalTable[i].DetectionCounter);
		xil_printf("SignalTable[%d] freq=%d power=%d len=%d per=%d delta_samples=%d delta_rf=%d\n\r", i, (u32)SignalTable[i].Freq, (s32)SignalTable[i].power, SignalTable[i].imp_len_max, SignalTable[i].imp_per_max, SignalTable[i].delta_samples, SignalTable[i].delta_rf);
	}


	return;

}
void SortByDecrease(u8 mode)
{
   int k;
   TSignal x;

   for (int i=0; i<Pelengator.SignalCnt-1; i++)
   {
        k = i;
        x =  (SignalTable[i]);

        for (int j=i+1; j<Pelengator.SignalCnt; j++)
        {
            //if (((sSignal*)SignalList->Items[j])->DetectionCounter > x->DetectionCounter)
            if (CompareSignals(&(SignalTable[j]), &x, mode))
            {
                k = j;
                x = (SignalTable[k]);
            };
        };
        SignalTable[k] = SignalTable[i];
        SignalTable[i] = x;

    };
   return;
};
//---------------------------------------------------------------------------
u8 CompareSignals(TSignal *y, TSignal *x, u8 method)
{
    // // сортировка с приоритетом по Position которые помечены !=-1 в начало списка
    // if( (y->Position != -1) && (x->Position != -1) )
    // {
    //     if(y->Position < x->Position)
    //         return 1;
    //     else
    //         return 0;
    // }
    // else if ( (y->Position != -1) && (x->Position == -1) )
    //     return 1;
    // else if ( (y->Position == -1) && (x->Position != -1) )
    //     return 0;
    // else
    {
        switch(method)
        {
        case 0:
        	//xil_printf("y->power=%d x->power=%d \n\r", y->power, x->power);
            // мощность
            if (y->power > x->power)
                return 1;
            else
                return 0;
            break;
        case 1:
            // счетчик
            if (y->DetectionCounter > x->DetectionCounter)
                return 1;
            else
                return 0;
            break;
        case 2:
            // Частота
            if (y->Freq > x->Freq)
                return 1;
            else
                return 0;
            break;
        // case 3:
        //     // Время
        //     if (CompareDateTime(y->DetectTime, x->DetectTime) == 1)
        //         return 1;
        //     else
        //         return 0;
        //     break;
        // case 4: break;
        default : break;
        }
    }
    return 0;
}
double ampl_spec_X7[250];
double vec_m[10];
double vec_D[10];

double FindSpectrumMax_PP(u8 AntNum)
{
	//u64 t1= TimemkS();
    int StartAnalyse = Pelengator.StartAnalyse ;// 0; // 0;
    int StopAnalyse = Pelengator.StopAnalyse; // round(Pelengator.Len*FreqBand/Fd); // SamplesNumber;///UT.iAntCount;

	// На входе имеем два массива комплексных чисел arr_cF7(опорный) и arr_cF9(сигнальный) длиной MLen = 256 отсчетов
	volatile comp *arr_cF7 = &fftpsd_rx2[StartAnalyse];
	volatile comp *arr_cF9 = &fftpsd_rx1[StartAnalyse];
	static u16 MLen = 250; // Pelengator.StopAnalyse - Pelengator.StartAnalyse;


	//double p = Power2log20(cabs(arr_cF7[135]));//-3
	//xil_printf("AntNum = %d p = %d freq = %d\r\n", AntNum, (s32)p, (u32)(lo_freq_hz/1e6));
	// 1. ВЫЧИСЛЕНИЕ ПОРОГА ОБНАРУЖЕНИЯ

	double Knorm = -130;
	double Pob_min = -60;

	//double ampl_spec_X7[MLen];

	int N = 10; // Количество частей, на которое разбивается спектр
	int n = MLen/N; // Количество отсчётов в одной части
	if( (global_MaxHoldRepeatCount==Pelengator.MaxHoldRepeatCount-1) && (global_AntCount==0) )
	{
		// Вычисление статистики
		//double vec_m[N];
		//double vec_D[N];


		for (int j = 0; j < N; ++j)
		{
			// Высичляем среднее значение
			vec_m[j] = 0;

			for (int k = 0; k < n; ++k)
			{
				if(cabs(arr_cF7[j*n + k]) != 0)
					ampl_spec_X7[j*n + k] = 20*log10(cabs(arr_cF7[j*n + k])) + Knorm;
				else
					ampl_spec_X7[j*n + k] = -240;
				vec_m[j] += ampl_spec_X7[j*n + k]/n;
			}

			// Вычисляем дисперсию
			vec_D[j] = 0;

			for (int k = 0; k < n; ++k)
			{
				vec_D[j] += pow(ampl_spec_X7[j*n + k] - vec_m[j], 2)/n;
			}
		}//for (int j = 0; j < N; ++j)

		double sensivity = Pelengator.sensivity;// 3;//5
		// Выбор наименьшего порога обнаружения
		Pob_min = vec_m[0] + sensivity*sqrt( vec_D[0] );

		for(int j = 1; j < N; ++j)
		{

			double Pob_j = vec_m[j] + sensivity*sqrt( vec_D[j] );

			if(Pob_j < Pob_min)
			{
				Pob_min = Pob_j;
			}
		}//for(int j = 1; j < N; ++j)
		Pob_min_array[global_ind]=Pob_min;
	}
	else
	{
		Pob_min = Pob_min_array[global_ind];
	}

	//xil_printf("Pob_min = %d \r\n", (s32)Pob_min);
	// 2. ОБНАРУЖЕНИЕ СИГНАЛОВ
	for(int i = 0; i < MLen; ++i)
	{
		if(ampl_spec_X7[i] > Pob_min)
		{
			int vec_indx[MLen];
			int count = 0;

			while(ampl_spec_X7[i] > Pob_min)
			{

				vec_indx[count] = i;
				count++;
				i++;

				if(i == MLen) break;
			}//while(ampl_spec_X7[i] > Pob_min)

			// Поиск максимума в полосе
			double P_max = -200;
			int indx_sign = -1; // Абсолютный индек положения максимума, он же индекс сигнала если выполнится проверка по полосе и по фазе внутри полосы
			int indx_sign_relative = -1; // Относительное значение индекса
			//double rf_max = 1000; // Значение разности фаз в максимуме
			comp rf_max = 0 + I*0;

			for(int j = 0; j < count; ++j)
			{
				if(ampl_spec_X7[vec_indx[j]] > P_max)
				{
					P_max = ampl_spec_X7[vec_indx[j]];
					indx_sign_relative = j;
					indx_sign = vec_indx[j];
					//rf_max = carg( arr_cF9[vec_indx[j]] / arr_cF7[vec_indx[j]] ) * 180.0/M_PI;
					rf_max = ( arr_cF9[vec_indx[j]] / arr_cF7[vec_indx[j]] );
				}
			}

			// Измерение полосы по уровню -10 дБ
				comp rf_j;
				double max_delta_rf = 0.0;
				int left_border = vec_indx[0]; // Левая граница полосы
				int right_border = vec_indx[count-1]; // Левая граница полосы
				double delta_rf = 0.0;

				// Движение влево от максимума
				for(int j = indx_sign_relative-1; j > -1; --j)
				{

					if( fabs(ampl_spec_X7[vec_indx[j]] - P_max) < 10 )
					{
						left_border = vec_indx[j];
						//double rf_j = carg( arr_cF9[vec_indx[j]] / arr_cF7[vec_indx[j]] ) * 180.0/M_PI;
						rf_j = ( arr_cF9[vec_indx[j]] / arr_cF7[vec_indx[j]] ) ;
						delta_rf = fabs(carg(rf_j / rf_max) * 180.0/M_PI);

						//if(fabs(rf_j - rf_max) > max_delta_rf)
						if( delta_rf > max_delta_rf )
						{
							max_delta_rf = delta_rf;// fabs(rf_j - rf_max);
						}
					}
					else
					{
						break;
					}
				}

				//double rf_j;

				// Движение вправо от максимума
				for(int j = indx_sign_relative+1; j < count; ++j)
				{

					if( fabs(ampl_spec_X7[vec_indx[j]] - P_max) < 10 )
					{
						right_border = vec_indx[j];
						//double rf_j = carg( arr_cF9[vec_indx[j]] / arr_cF7[vec_indx[j]] ) * 180.0/M_PI;
						rf_j = ( arr_cF9[vec_indx[j]] / arr_cF7[vec_indx[j]] );
						delta_rf = fabs(carg(rf_j / rf_max) * 180.0/M_PI);

						//if(fabs(rf_j - rf_max) > max_delta_rf)
						if( delta_rf > max_delta_rf )
						{
							max_delta_rf = delta_rf;//fabs(rf_j - rf_max);
						}
					}
					else
					{
						break;
					}
				}
				//if(max_delta_rf > 180) max_delta_rf-=180;

			// Проверка на отклонение разности фаз от нормы
				if( ((right_border - left_border + 1) >= Pelengator.max_delta_samples) && (max_delta_rf < Pelengator.max_delta_rf) )  // 1  50.0
				{
					// Добавляем сигнал в таблицу, индекс обнаруженного сигнала indx_sign
					/*P_max мощность
					indx_sign частота
					arr_cF7(опорный) и arr_cF9(сигнальный)
					left_border // Левая граница полосы
					right_border
					Pob_min// порог
					*/
					//xil_printf("indx_sign=%d\r\n", indx_sign);
					//xil_printf("max_delta_samples=%d max_delta_rf=%d ", right_border - left_border + 1, (s32)max_delta_rf);
					AddSignal_PP(                                            //
							FFTcodeToFreq_1(indx_sign+StartAnalyse)+lo_freq_hz/1e6,       //  double Freq,
							(P_max),                              //  double Power,
							arr_cF9[indx_sign],                              //  comp CompCh1,
							arr_cF7[indx_sign],                              //  comp CompCh2,
							AntNum,                                          //  int AntNum,
							FFTcodeToFreq_1(right_border+StartAnalyse)+lo_freq_hz/1e6,    //  double freq_max,
							FFTcodeToFreq_1(left_border+StartAnalyse)+lo_freq_hz/1e6,     //  double freq_min)
							(u16)(right_border - left_border + 1),
							(s8)max_delta_rf
							);


				}
		}//if(ampl_spec_X7[i] > Pob_min)
	}//for(int i = 0; i < MLen; ++i)

	//xil_printf("time find max = %d \r\n", (u32)(TimemkS() - t1));
	return Pob_min;

}
double Power2log20(double power)
{
	if(power != 0)
		return (20*log10(power)-100-30);
	else
		return -240;
}
int AddSignal_PP(
		double Freq,
		double Power,
		comp CompCh1,
		comp CompCh2,
		int AntNum,
		double freq_max,
		double freq_min,
		u16 delta_samples,
		s8 delta_rf
		)
{
//xil_printf("AddSignal_PP AntNum=%d Freq=%d Power=%d", AntNum, (s32)Freq,(s32)Power);
	//if (CheckSpurs(Freq,Power)>0)
	{
		struct TSignal Signal;
		memset(&Signal, 0, sizeof(TSignal));
      //Signal = new sSignal;
      Signal.Freq = Freq;
      //if(AntNum==0)
      {
    	  Signal.Power = Power;
    	  Signal.power = (s8)round(Power);
    	  Signal.fc_max = (u16)round(freq_max);
    	  Signal.fc_min = (u16)round(freq_min);
      }
      Signal.Sch1[AntNum] = CompCh1;
      Signal.Sch2[AntNum] = CompCh2;
      Signal.AntNum[AntNum] = AntNum;
      Signal.PelengFi = 0;
      Signal.PelengTeta = 0;
      Signal.delta_samples = delta_samples;
      Signal.delta_rf = delta_rf;

      Signal.imp_len_max = (u32)(Samples2Time_mks(dsp_out_params_arr[global_ind].duration)*100.0);
      Signal.imp_len_min = (u32)(Samples2Time_mks(dsp_out_params_arr[global_ind].duration)*100.0);
      Signal.imp_per_max = (u32)(Samples2Time_mks(dsp_out_params_arr[global_ind].period)*100.0);
      Signal.imp_per_min = (u32)(Samples2Time_mks(dsp_out_params_arr[global_ind].period)*100.0);




      //xil_printf("Signal.power = %d \r\n", Signal.power);
      volatile u16 i=0;

      for (i=0; i< Pelengator.SignalCnt;i++)
      {
    	  TSignal * TmpSignal=&(SignalTable[i]);

    	bool btmp = comparison_signal_fz_band(
    		  TmpSignal->Freq, //  double fz_1,
			  Freq, // double fz_2,
			  fabs(Signal.fc_max - Signal.fc_min), // double band_1,
			  fabs(freq_max - freq_min), // double band_2,
			  Pelengator.porog_fz,
			  Pelengator.porog_band);
    	  if (btmp == true)
    	  //float Tmp = fabs(Freq - TmpSignal->Freq);

    	 // if (Tmp < 3)
    	  {
         	TmpSignal->DetectTime = TimeS();
         	//if(AntNum==0)
         	{
         		TmpSignal->Power = Power;
            	TmpSignal->power = (s8)round(Power);
            	TmpSignal->fc_max = (u16)round(freq_max);
            	TmpSignal->fc_min = (u16)round(freq_min);
         	}
            //TmpSignal->imp_len_max = TmpSignal->DetectionCounter;
            TmpSignal->DetectionCounter++;
            TmpSignal->Sch1[AntNum] = Signal.Sch1[AntNum];
            TmpSignal->Sch2[AntNum] = Signal.Sch2[AntNum];
            TmpSignal->AntNum[AntNum] = AntNum;
            TmpSignal->AntNumCnt++;
            TmpSignal->Refresh=1;
            TmpSignal->delta_samples = delta_samples;
            TmpSignal->delta_rf = delta_rf;

            TmpSignal->imp_len_max = (u32)(Samples2Time_mks(dsp_out_params_arr[global_ind].duration)*100.0);
            TmpSignal->imp_len_min = (u32)(Samples2Time_mks(dsp_out_params_arr[global_ind].duration)*100.0);
            TmpSignal->imp_per_max = (u32)(Samples2Time_mks(dsp_out_params_arr[global_ind].period)*100.0);
            TmpSignal->imp_per_min = (u32)(Samples2Time_mks(dsp_out_params_arr[global_ind].period)*100.0);



           //xil_printf("find!\r\n");

            return 1;
    	  }
      }

      Signal.DetectTime = TimeS();
      Signal.DetectionCounter=1;
      //Signal.imp_len_max = 1;

      AddToSignalTable(Signal);
     // xil_printf("new! %d\r\n", Pelengator.SignalCnt);
      return 1;
   }
   //else {return 0;}
}
bool comparison_signal_fz_band(double fz_1, double fz_2, double band_1, double band_2, double porog_fz, double porog_band)
{
    // porog_fz = 1 МГц (можно регулировать) - допустимая разность частот
	// porog_band = 0.6 - допустимое отношение полос (60 %)

	if(fabs(fz_1 - fz_2) > porog_fz) {
        return false;
    }

	if(band_1/band_2 > 1.0)
	{
		if(band_2/band_1 < porog_band)
		{
			return false;
		}
	}
	else
	{
		if(band_1/band_2 < porog_band)
		{
			return false;
		}
	}

	return true;
}
u32 FindMax_Scan(u8* buf)
{

	volatile u16 CmdNum = *(u16*)(buf + PAYLOAD_CMD_OFFSET);
	volatile u8 *DataBuf = (buf + 20);
	volatile u32 splencode = (*(u16 *)(DataBuf + 0));
	splencode &= ~(0x3);
	splencode >>= 1;

	volatile u8 aFDcount = DataBuf[4];

	volatile u16 aTimerPeriod = (*(u16 *)(DataBuf + 5));
	volatile u16 aDelayCycle = (*(u16 *)(DataBuf + 7));

	volatile u16 aDetectThreshold = (*(u16 *)(DataBuf + 9));
	Pelengator.Threshold = aDetectThreshold;

	volatile u8 aiAntCount = (*(u8 *)(DataBuf + 11 ));
	volatile u8 aiTaskNumber =  (*(u8 *)(DataBuf + 12));
	volatile u8 aiExternalSyncEn = (*(u8 *)(DataBuf + 13));
	volatile u8 pulse_meas_en = 1;


	volatile u16 i = 0;
	volatile u16 offset = 14;
	for(i = 0;i<aFDcount;i++)
	{
	   Pelengator.aiFmin[i] = (*(u16 *)(DataBuf + offset+0 + i*5)); // DataBuf[offset+1 + i*5] | (DataBuf[offset + i*4] << 8);
	   Pelengator.aiFCount[i] = (*(u16 *)(DataBuf + offset+2 + i*5)); // DataBuf[offset+3 + i*5] | (DataBuf[offset+2 + i*4] << 8);
	   Pelengator.RxGain_Array[i] = (*(u8 *)(DataBuf + offset+4 + i*5)); //
	}

	volatile u8 TaskRepeatCount = (*(u8 *)(DataBuf + offset + aFDcount*(5)));
	volatile u8 FreqRepeatCount = (*(u8 *)(DataBuf + offset + aFDcount*(5) + 1));
	volatile u8 MaxHoldRepeatCount = (*(u8 *)(DataBuf + offset + aFDcount*(5) + 2));
	volatile u8 locMaxHoldRepeatCount = MaxHoldRepeatCount;
	volatile u16 MaxHoldDelay = (*(u16 *)(DataBuf + offset + aFDcount*(5) + 3));
	volatile u16 DetectFreq = (*(u16 *)(DataBuf + offset + aFDcount*(5) + 5));
	volatile u16 DetectFreqDelta = (*(u16 *)(DataBuf + offset + aFDcount*(5) + 7));

	Pelengator.DetectFreq = DetectFreq;
	Pelengator.DetectFreqDelta = DetectFreqDelta;

	volatile u8 locFreqRepeatCount = FreqRepeatCount;
	if(FreqRepeatCount == 0)
	{
	   FreqRepeatCount = 1;
	}

	Pelengator.Len = splencode;
	////--------------
	CmdNum = ETH_CMD_WB_SCAN;
	//volatile u32 FreqStepNumber = Pelengator.FreqCount; //  116;
	volatile u32 DataSize = 0; // (Pelengator.Len << 3) * FreqStepNumber ;
	volatile u32 BytesSentOverall = 0;
	volatile int idx;




	volatile u8 FDcountINDEX = 0, ANTcountINDEX = 0, DetectEn = 0;
	volatile u32 BytesSentReal = 0, DetectRezultFirstAnt = 0;
	volatile s32 DetectRezult = 0;
	volatile u16 AttenFreqIndex = 0;
	volatile u16 AntAttenFreqIndex = 0;
	volatile u8 local_rx_gain = 0;


	DataSize = 0;
	u16 NStep = 0;

	for(FDcountINDEX=0;FDcountINDEX < aFDcount; FDcountINDEX++)
	{
		DataSize += ((Pelengator.Len << Pelengator.ch_num_shift)+Pelengator.extra_byte_fft) * Pelengator.aiFCount[FDcountINDEX] * FreqRepeatCount  * MaxHoldRepeatCount * aiAntCount;
		NStep = Pelengator.aiFCount[FDcountINDEX];
	}

	Xil_DCacheInvalidateRange((UINTPTR)(RX_BUFFER_BASE),  DataSize);
	while(TaskRepeatCount--)
	{



		BytesSentReal = 0;
		BytesSentOverall = 0;
		AttenFreqIndex = 0;
		AntAttenFreqIndex = 0;
		global_ind = AttenFreqIndex;



		for(FDcountINDEX=0;FDcountINDEX < aFDcount; FDcountINDEX++)
		{


			for(i=0;i < Pelengator.aiFCount[FDcountINDEX]; i++)
			{
				lo_freq_hz = lo_freq_hz_next;

				if( i == (Pelengator.aiFCount[FDcountINDEX]-1) )
				{
					if( FDcountINDEX == (aFDcount-1) )
					{
						lo_freq_hz_next = Pelengator.aiFmin[0]*1e6;
						// xil_printf("FDcountINDEX == (aFDcount-1) \r\n");
					}
					else
					{
						lo_freq_hz_next = Pelengator.aiFmin[FDcountINDEX+1]*1e6;
						//   xil_printf("FDcountINDEX != (aFDcount-1) \r\n");
					}
				}
				else
				{
				   lo_freq_hz_next = lo_freq_hz+Pelengator.CalibrationFreqStep;
				   //   xil_printf("i != (Pelengator.aiFCount[FDcountINDEX]-1) \r\n");
				}

				locFreqRepeatCount = FreqRepeatCount;

				while(locFreqRepeatCount--)
				{
					global_FreqRepeatCount = locFreqRepeatCount;

					locMaxHoldRepeatCount = MaxHoldRepeatCount;

					while(locMaxHoldRepeatCount--)
					{
						global_MaxHoldRepeatCount = locMaxHoldRepeatCount;

						for(ANTcountINDEX=0;ANTcountINDEX < aiAntCount; ANTcountINDEX++)
						{
							global_AntCount=ANTcountINDEX;
							RePack((UINTPTR)(RX_BUFFER_BASE+BytesSentOverall), AttenFreqIndex, AntAttenFreqIndex);//AntAttenFreqIndex
							FindSpectrumMax_PP(ANTcountINDEX);
							BytesSentOverall = BytesSentOverall + (Pelengator.Len << Pelengator.ch_num_shift);// + Pelengator.extra_byte_fft;
							AntAttenFreqIndex++;
						}//for(ANTcountINDEX=0;ANTcountINDEX < aiAntCount; ANTcountINDEX++)


					}// while(locMaxHoldRepeatCount--)

				}//while(locFreqRepeatCount--)

			   AttenFreqIndex++;
			   global_ind = AttenFreqIndex;

			}//for(i=0;i < Pelengator.aiFCount[FDcountINDEX]; i++)

		}//for(FDcountINDEX=0;FDcountINDEX < aFDcount; FDcountINDEX++)

	}//while(TaskRepeatCount--)


	return 0;
}
