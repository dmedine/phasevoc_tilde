

#include "m_pd.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

static t_class *phasevoc_tilde_class;

typedef struct _phasevoc_tilde{
  t_object x_obj;
  t_sample f;
  t_sample *circBuff, *sampHoldBuff;
  long sampleRate, n, buffer_limit;
  long circBuffIndex, delayCount;
  long sampHoldBuffLength, sampHoldBuffIndex, circBuffToSampHoldCopyIndex;
  long bufferPosition;
  t_float circBuffLength;
  t_float sampHoldBuffWindowIndex, sampHoldBuffWindowStart;
  t_float circBuffWindowFrontIndex, circBuffWindowBackIndex, circBuffWindowStart;
  t_float sampHoldBuffWindowFrontIndex, sampHoldBuffWindowBackIndex;
  t_float timestretch, pitchshift;
  t_float circBuffSeconds;
  t_float twoPi;
  long hopSize, fftSize, fftHalfSize, dsp_tick;
  long inputTime, outputTime, overlap, overlap_chunk;
  t_sample *inFrontBuff, *inBackBuff, *outBuff;
  t_sample *inFrontShift, *inBackShift;
  t_sample *nonoverlappedOutBuff;
  t_sample *analWindow, *synthWindow;
  t_sample *frontReal,*frontImag, *backReal, *backImag;
  t_sample *prevReal, *prevImag, *conjReal,  *conjImag;
  t_sample *outReal, *outImag, *outSpectra, *rsqrt;
  t_int fft_tick, overlapChunk;
  t_int phaselock, pause, waszeroflag;

  t_outlet *buffPos;

  
 

}t_phasevoc_tilde;


/***********function declarations*******************/
float phasevoc_tilde_interpolate(t_float *buffer, long bufferLength, double findex);
t_float phasevoc_tilde_cents(t_phasevoc_tilde *x);
void phasevoc_tilde_samphold(t_phasevoc_tilde *x, t_floatarg g);
void phasevoc_tilde_phaselock(t_phasevoc_tilde *x, t_floatarg g);
void phasevoc_tilde_hopsize(t_phasevoc_tilde *x, t_floatarg g);
void phasevoc_tilde_windowsize(t_phasevoc_tilde *x, t_floatarg g);
void phasevoc_tilde_buffsize(t_phasevoc_tilde *x, t_floatarg g);
void phasevoc_tilde_overlap(t_phasevoc_tilde *x, t_floatarg g);
void phasevoc_tilde_reinit(t_phasevoc_tilde *x);
void phasevoc_tilde_pause(t_phasevoc_tilde *x);
void phasevoc_tilde_unpause(t_phasevoc_tilde *x);
static void init_window(t_phasevoc_tilde *x);



/************utility functions***************/

/********spline interpolation*********/
float phasevoc_tilde_interpolate(t_float *buffer, long bufferLength, double findex)
{
  long lindex = findex;
  t_sample fr = findex - lindex;

  float p0=buffer[(lindex-2) % bufferLength];
  float p1=buffer[(lindex-1) % bufferLength];
  float p2=buffer[(lindex) % bufferLength];
  float p3=buffer[(lindex+1) % bufferLength];
  float p4=buffer[(lindex+2) % bufferLength];
  float p5=buffer[(lindex+3) % bufferLength];

  return p2 + 0.04166666666*fr*((p3-p1)*16.0+(p0-p4)*2.0
				+ fr *((p3+p1)*16.0-p0-p2*30.0- p4
				       + fr *(p3*66.0-p2*70.0-p4*33.0+p1*39.0+ p5*7.0- p0*9.0
					      + fr *( p2*126.0-p3*124.0+p4*61.0-p1*64.0- p5*12.0+p0*13.0
						      + fr *((p3-p2)*50.0+(p1-p4)*25.0+(p5-p0)*5.0)))));
}

/***************cents*******************/
t_float phasevoc_tilde_cents(t_phasevoc_tilde *x)
{
  t_float g;
  g = (x->pitchshift * .01) + 69;
  g = 8.17579891564 * exp(.0577622650 * g);
  g = g/440;
  return g;
}

/*************process replacing************/
t_int *phasevoc_tilde_perform (t_int *w)
{
  t_phasevoc_tilde *x = (t_phasevoc_tilde *)(w[1]);
  t_sample  *in = (t_sample *)(w[2]);
  t_sample *out1 = (t_sample *)(w[3]);
  int n = (int)(w[4]);
  t_float n2 = (t_float)n;
  x->n = n;
  int i, j;
  t_float psfact = phasevoc_tilde_cents(x);
  t_float tsfact = x->timestretch;
  t_float buffpos;
  long maskFFT=x->fftSize-1;
 

  //unless samphold or pause is on, always fill the circBuff
  if(x->pause!=1)
    {
      for(i=0; i<n; i++)
	{
	  x->circBuff[x->circBuffIndex++] = in[i];
	  while(x->circBuffIndex > (x->circBuffLength-1))x->circBuffIndex-=(x->circBuffLength);
	  while(x->circBuffIndex < 0 )x->circBuffIndex += x->circBuffLength;
	  
	}
      x->delayCount+=n;
    }

 
 if(x->delayCount>x->fftSize)
    {
      if(!x->sampHoldBuffLength)
	{
	  //init the indeces for this read-out fromt the circ buff
	  x->circBuffWindowBackIndex = x->circBuffWindowStart;
	  x->circBuffWindowFrontIndex = x->circBuffWindowBackIndex + (x->hopSize*psfact);
	}

      else
	{
	  //init the indeces for this read-out fromt the sampholdbuff
	  x->sampHoldBuffWindowBackIndex = x->sampHoldBuffWindowStart;
	  x->sampHoldBuffWindowFrontIndex = x->sampHoldBuffWindowBackIndex + (x->hopSize*psfact);
	}


      //check for dsp_tick
      if(x->dsp_tick == (x->buffer_limit*x->overlap))
	x->dsp_tick = 0;
      
      
      if((x->dsp_tick%x->buffer_limit)==0)
	{
	  //buffer in the samples if no sampbuff
	  if(!x->sampHoldBuffLength)
	    {
	      for(i=0; i<x->fftSize; i++)
		{
		  x->inFrontBuff[i] =  phasevoc_tilde_interpolate(x->circBuff, x->circBuffLength, x->circBuffWindowFrontIndex);
		  x->inBackBuff[i] =  phasevoc_tilde_interpolate(x->circBuff, x->circBuffLength, x->circBuffWindowBackIndex);
		  
		  //increment the indices by the pitchshifting factor
		  x->circBuffWindowFrontIndex+=psfact;
		  x->circBuffWindowBackIndex+=psfact;
		  
		  //wrap 'em
		  while(x->circBuffWindowFrontIndex>(x->circBuffLength-1))
		      x->circBuffWindowFrontIndex-=(x->circBuffLength);
	
		  while(x->circBuffWindowFrontIndex<0)
		      x->circBuffWindowFrontIndex+=x->circBuffLength;
		  while(x->circBuffWindowBackIndex>(x->circBuffLength-1))x->circBuffWindowBackIndex-=(x->circBuffLength);
		  while(x->circBuffWindowBackIndex<0)x->circBuffWindowBackIndex+=x->circBuffLength;
		}
	    }
	  //buffer in the samples if sampbuff
	  else
	    {
	      for(i=0; i<x->fftSize; i++)
		{
		  x->inFrontBuff[i] =  phasevoc_tilde_interpolate(x->sampHoldBuff, x->sampHoldBuffLength, x->sampHoldBuffWindowFrontIndex);
		  x->inBackBuff[i] =  phasevoc_tilde_interpolate(x->sampHoldBuff, x->sampHoldBuffLength, x->sampHoldBuffWindowBackIndex);
		  
		  //increment the indices by the pitchshifting factor
		  x->sampHoldBuffWindowFrontIndex+=psfact;
		  x->sampHoldBuffWindowBackIndex+=psfact;
		  
		  //wrap 'em
		  while(x->sampHoldBuffWindowFrontIndex>(x->sampHoldBuffLength-1))
		      x->sampHoldBuffWindowFrontIndex-=(x->sampHoldBuffLength);
		   
		  while(x->sampHoldBuffWindowFrontIndex<0)
		      x->sampHoldBuffWindowFrontIndex+=x->sampHoldBuffLength;
		    
		  while(x->sampHoldBuffWindowBackIndex>(x->sampHoldBuffLength-1))x->sampHoldBuffWindowBackIndex-=(x->sampHoldBuffLength);
		  while(x->sampHoldBuffWindowBackIndex<0)x->sampHoldBuffWindowBackIndex+=x->sampHoldBuffLength;
		}
	    }

	  //window the signal
	  for(i=0; i<x->fftSize; i++)
	    {

	      x->inFrontBuff[i] *= x->analWindow[i];
	      x->inBackBuff[i] *= x->analWindow[i];

	    }
	  
	  
	  
	  //go to frequency domain
	  mayer_realfft(x->fftSize, x->inFrontBuff);
	  mayer_realfft(x->fftSize, x->inBackBuff);

	  
	  
	  //unpack the bizzarely packed mayer ffts
	  for(i=0; i<=x->fftHalfSize; i++)// halfSizeFFT = nyquist
	    {
	      x->frontReal[i] = x->inFrontBuff[i];
	      x->backReal[i] = x->inBackBuff[i];

	    }
	  x->frontImag[0] = x->backImag[0]= 0;  // 0 DC
	  
	  for(i=(x->fftSize-1), j=1; i>x->fftHalfSize; i--, j++)
	    {

	      x->frontImag[j] = x->inFrontBuff[i];
	      x->backImag[j] = x->inBackBuff[i];


	    }
	  
	  x->frontImag[x->fftHalfSize]=x->backImag[x->fftHalfSize]=0; //halfSizeFFT empty
	  
	  for(i = 0; i <= x->fftHalfSize; i++)
	    {
	      //first normalize over magnitude so we get phase only
	      x->rsqrt[i] = q8_rsqrt((x->prevReal[i] * x->prevReal[i]) +
				     (x->prevImag[i] * x->prevImag[i]) +
				     .000000000000000000001f);
	      x->prevReal[i] *= x->rsqrt[i];
	      x->prevImag[i] *= x->rsqrt[i];

	      x->conjReal[i] = (x->prevReal[i] * x->backReal[i]) + (x->prevImag[i] * x->backImag[i]) +
	      	.0000000000000001;
	      x->conjImag[i] = (x->prevImag[i] * x->backReal[i]) - (x->prevReal[i] * x->backImag[i]);

	    	      
	      x->rsqrt[i] = q8_rsqrt((x->conjReal[i] * x->conjReal[i]) +
	      			     (x->conjImag[i] * x->conjImag[i]));
	      x->conjReal[i] *= x->rsqrt[i];
	      x->conjImag[i] *= x->rsqrt[i];

	      //then add the front window (complex multiply) this time keeping magnitudes, and store the info for the next time around
	      x->outReal[i] = (x->conjReal[i] * x->frontReal[i]) - (x->conjImag[i] * x->frontImag[i]);
	      
	      x->outImag[i] = (x->conjReal[i] * x->frontImag[i]) + (x->conjImag[i] * x->frontReal[i]);
	      x->prevImag[i] = x->outImag[i];
	      x->prevReal[i] = x->outReal[i];
	    }

	  //phase has been propogated, repack for ifft
	  for(i=0; i<=x->fftHalfSize; i++)  // +1 to include Nyquist
	    x->outSpectra[i] = x->outReal[i];
	  
	  for(j = x->fftHalfSize -1, i = x->fftHalfSize + 1; i < x->fftSize; j--, i++)
	    x->outSpectra[i] = x->outImag[j];
	  
	  //back to time domain
	  mayer_realifft(x->fftSize, x->outSpectra);

	  //william brent's overlap/add function:
	  //window the synthesis
	  for(i=0; i<x->fftSize; i++)
	    x->outSpectra[i] *= x->synthWindow[i];//synthWindow is 2/3/fftSize*analWindow
	  
	  //  first shift output in nonoverlapped buffer
	  for(i=0; i<((x->overlap-1)*x->fftSize); i++)
	    x->nonoverlappedOutBuff[i] = x->nonoverlappedOutBuff[x->fftSize+i];
	 
	  //  then, write a new window in
	  for(i=0; i<x->fftSize; i++)
	    x->nonoverlappedOutBuff[((x->overlap-1)*x->fftSize)+i] = x->outSpectra[i];
	  
	  // init this chunk of the final output so it can be summed in the for() below
	  for(i=0; i<(x->hopSize); i++)
	    x->outBuff[i] = 0.0;
	  
	  // do the overlap/add: add the last hopsize of the first
	  for(i=0; i<x->overlap; i++)
	    for(j=0; j<(x->hopSize); j++)
	      x->outBuff[j] += x->nonoverlappedOutBuff[(i*x->fftSize)+((x->overlap-i-1)*x->hopSize)+j/*overlap_chunk*/];
	  x->fft_tick=0;
	}
      
      //output
      for(i=0; i<n; i++, out1++)
	*out1 = x->outBuff[(x->fft_tick*n)+i];
      if(!x->sampHoldBuffLength)
	buffpos = 1.0-(x->circBuffIndex - x->circBuffWindowStart) / x->circBuffLength;
      else
	buffpos = 1.0-(x->sampHoldBuffIndex - x->sampHoldBuffWindowStart) / x->sampHoldBuffLength;
      while(buffpos>=1.0)buffpos -= 1.0;
      while(buffpos<=0.0)buffpos += 1.0;
      outlet_float(x->buffPos, buffpos);
      x->fft_tick++;
      x->dsp_tick++;

      //increment the read start point and wrap it
      if(!x->sampHoldBuffLength)
	{
	  x->circBuffWindowStart+=n*tsfact;
	  while(x->circBuffWindowStart > (x->circBuffLength-1))x->circBuffWindowStart-=x->circBuffLength;
	  while(x->circBuffWindowStart<0)x->circBuffWindowStart+=x->circBuffLength;
	}

      else
	{
	  x->sampHoldBuffWindowStart+=n*tsfact;
	  while(x->sampHoldBuffWindowStart > (x->sampHoldBuffLength-1))x->sampHoldBuffWindowStart-=x->sampHoldBuffLength;
	  while(x->sampHoldBuffWindowStart<0)x->sampHoldBuffWindowStart+=x->sampHoldBuffLength;
	}

      return(w+5);
    }
 else
   {
     *out1++ = 0.0;
     return(w+5);

   }
   
}

//inlet functions-------------------------------------------------------------------
/***********timestretch***********/
void phasevoc_tilde_timestretch (t_phasevoc_tilde *x, t_floatarg g)
{
 
  x->timestretch = g;

}


/***********pitchshift***********/
void phasevoc_tilde_pitchshift (t_phasevoc_tilde *x, t_floatarg g)
{

   x->pitchshift = g;
}


//creation argument functions--------------------------------------------------------
/************************windowsize*****************/
void phasevoc_tilde_windowsize (t_phasevoc_tilde *x, t_floatarg g)
{
   float base;
   float y = g;
  
    y = log(y);
    base = log(2);
    y = y/base;
   int inty = (int)y;
    if(y==inty)
      {
	if(g>=256)
	  {
	    x->fftSize = g;
	    post("windowsize = %f", g);
	  }
	else
	  error("windowsize must be >= 256: windowsize = 1024");
      }

    else
      {
	error("windowsize must be a power of 2: windowsize = 1024");
	x->fftSize = 1024;
      }
 
}

/***********************hopsize******************/
void phasevoc_tilde_hopsize (t_phasevoc_tilde *x, t_floatarg g)
{

if(x->fftSize % (int)g != 0)
    {
      error("windowsize must be divisible by hopsize: hopsize = %f", (float)x->fftSize/4);
        x->hopSize = x->fftSize/4;
    }

 else if(g<64)
   {
     error("minimum hopsize right now is 64, sorry");
     x->hopSize = 256;
   }

else
  {
    x->hopSize = g;
    post("hopsize = %d", g);
  }

}

/***********************buffsize**********************/
void phasevoc_tilde_buffsize (t_phasevoc_tilde *x, t_floatarg g)
{
  if (g<=0)
    {
      post("buffsize must be greater than 0: buffsize = 5");
      x->circBuffSeconds = 5;
    }
  else
    {
      x->circBuffSeconds = g;
      // memset(x->circBuff, 0, sizeof(float) * x->circBuffLength);
      x->circBuffLength = x->circBuffSeconds * x->sampleRate;
      
      x->circBuff = (t_float *)t_getbytes(0);
      x->circBuff = (t_float *)t_resizebytes(x->circBuff, 0, sizeof(float) * x->circBuffLength);
      
      phasevoc_tilde_reinit(x);
      
      post("input buffer is %f seconds", g);
    }
}

/**************overlap******************/
void phasevoc_tilde_overlap(t_phasevoc_tilde *x, t_floatarg g)
{}

//functions by messages------------------------------------------


/*****************samphold****************/
void phasevoc_tilde_samphold (t_phasevoc_tilde *x, t_floatarg g)
{
  long dummy, dummy2, dummy3;
  int i;
  dummy = x->circBuffIndex;
  if (g<=0.0)
    {
      x->sampHoldBuffLength = 0;
      x->sampHoldBuffIndex = 0;
      x->sampHoldBuffWindowIndex = 0;
      x->sampHoldBuffWindowStart = 0;
      phasevoc_tilde_reinit(x);
    }
  else
    {
      x->pause = 1;
      
      x->sampHoldBuffLength = x->sampleRate * g;
      if(x->sampHoldBuffLength>=x->circBuffLength)error("sh interval must be less then %f seconds", x->circBuffSeconds);
      else post("sample and hold interval is %f seconds", g);
      //copy to sampbuff the last g seconds of circbuff
    
      x->sampHoldBuffWindowFrontIndex = x->sampHoldBuffWindowBackIndex = x->sampHoldBuffWindowStart = x->sampHoldBuffIndex = 0;

     
      if(x->circBuffIndex - x->sampHoldBuffLength>=0)
	{
	  dummy = x->circBuffIndex-x->sampHoldBuffLength;
	  memcpy(x->sampHoldBuff, x->circBuff + (dummy), sizeof(float)*x->sampHoldBuffLength);
	}
      else
	{
	  dummy = x->circBuffIndex - 0;
	  dummy2 = x->sampHoldBuffLength-dummy;
	  dummy3 = x->circBuffLength - dummy2;

	  memcpy(x->sampHoldBuff, x->circBuff + dummy3, sizeof(float)*dummy2);
	  memcpy(x->sampHoldBuff + dummy2, x->circBuff, sizeof(float)*dummy);
	}
    }

}

/**********************phaselock**********************/
void phasevoc_tilde_phaselock(t_phasevoc_tilde *x, t_floatarg g)
{
  int i;

  if (g != 0)
    x->phaselock = 1;
    
  else
    x->phaselock = 0;
}

/**************reinit*****************/
void phasevoc_tilde_reinit(t_phasevoc_tilde *x)
{
  int i;
  for(i=0; i<x->circBuffLength; i++)
    x->circBuff[i] = 0.0;
  //memset(x->circBuff, 0, sizeof(float) * x->circBuffLength);
  x->circBuffIndex = 0;
  x->circBuffWindowBackIndex = x->circBuffWindowStart = 0;
  x->sampHoldBuffIndex = x->sampHoldBuffWindowIndex = x->sampHoldBuffWindowStart = x->circBuffToSampHoldCopyIndex = 0;
  x->sampHoldBuffLength = 0;
  phasevoc_tilde_pause(x);

  post("phasevoc~ reinitialized\n");

}

/****************pause***************/
void phasevoc_tilde_pause(t_phasevoc_tilde *x)
{
  x->pause = 1;
  post("paused");
}

/****************unpause***************/
void phasevoc_tilde_unpause(t_phasevoc_tilde *x)
{
  x->pause = 0;
  post("unpaused");
}


//------------------------------------------------------------------------
//------------------------------------------------------------------------
//------------------------------------------------------------------------
/************new***********/
void *phasevoc_tilde_new(t_symbol *s, int argc, t_atom *argv)
{
 

  t_phasevoc_tilde *x = (t_phasevoc_tilde *)pd_new(phasevoc_tilde_class);
 
 

  inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("timestretch"));
  inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("pitchshift"));
  outlet_new(&x->x_obj, &s_signal);
  x->buffPos = outlet_new(&x->x_obj, &s_float);

  x->sampleRate = sys_getsr();
  x->circBuffSeconds = 5;
  x->circBuffLength = x->sampleRate * x->circBuffSeconds;
  x->circBuffIndex = 0;
  x->circBuffWindowBackIndex = x->circBuffWindowStart = 0;
  x->sampHoldBuffIndex = x->sampHoldBuffWindowIndex = x->sampHoldBuffWindowStart = x->circBuffToSampHoldCopyIndex = 0;
  x->sampHoldBuffLength = 0;
  x->timestretch = 1.0;
  x->pitchshift = 0.0;
  x->fftSize = 1024;
  x->fftHalfSize = x->fftSize >> 1;
  x->hopSize = 64;
  x->overlapChunk = x->fftSize>>2;
  x->overlap = x->fftSize/x->hopSize;
  x->twoPi = 8.0f * atanf(1.0f);
  x->bufferPosition = 0;
  x->n = 64;
  x->inputTime = x->outputTime = x->dsp_tick = 0;
  x->buffer_limit = x->fftSize/x->n/x->overlap;
  x->fft_tick = 0;
  x->outputTime = 0;
  x->phaselock = 0;
  x->pause = 0;
  x->waszeroflag = 0;



  while (argc > 0)
    {
      t_symbol *firstarg = atom_getsymbolarg(0, argc, argv);
      if (!strcmp(firstarg->s_name, "-buffsize"))
        {
	  phasevoc_tilde_buffsize(x, atom_getfloatarg(1, argc, argv));
	  argc-=2, argv+=2;
        }
      else if (!strcmp(firstarg->s_name, "-windowsize"))
        {
	  phasevoc_tilde_windowsize(x, atom_getfloatarg(1, argc, argv));
	  argc-=2, argv+=2;
        }
      else if (!strcmp(firstarg->s_name, "-hopsize"))
        {
	  phasevoc_tilde_hopsize(x, atom_getfloatarg(1, argc, argv));
	  argc-=2, argv+=2;
        }
      else
	{
	  pd_error(x, "phasevoc: %s: unknown flag or argument missing",
		   firstarg->s_name);
	  argc--, argv++;
	}

    }


  
  //
  x->circBuff = (t_float *)t_getbytes(0);
  x->circBuff = (t_float *)t_resizebytes(x->circBuff, 0, sizeof(float) * x->circBuffLength);
  //
  x->sampHoldBuff = (t_float *)t_getbytes(0);
  x->sampHoldBuff = (t_float *)t_resizebytes(x->sampHoldBuff, 0, sizeof(float) * x->circBuffLength);
  //
  x->inFrontBuff = (t_float *)t_getbytes(0);
  x->inFrontBuff = (t_float *)t_resizebytes(x->inFrontBuff, 0, sizeof(float) * x->fftSize);
  //
  x->inBackBuff = (t_float *)t_getbytes(0);
  x->inBackBuff = (t_float *)t_resizebytes(x->inBackBuff, 0, sizeof(float) * x->fftSize);
  //
  x->outBuff = (t_float *)t_getbytes(0);
  x->outBuff = (t_float *)t_resizebytes(x->outBuff, 0, sizeof(float) * (x->hopSize + x->n));
  //
  x->inFrontShift = (t_float *)t_getbytes(0);
  x->inFrontShift = (t_float *)t_resizebytes(x->inFrontShift, 0, sizeof(float) * x->fftSize);
  //
  x->inBackShift = (t_float *)t_getbytes(0);
  x->inBackShift = (t_float *)t_resizebytes(x->inBackShift, 0, sizeof(float) * x->fftSize);
  //
  x->frontReal = (t_float *)t_getbytes(0);
  x->frontReal = (t_float *)t_resizebytes(x->frontReal, 0, sizeof(float) * (x->fftHalfSize+1));
  //
  x->frontImag = (t_float *)t_getbytes(0);
  x->frontImag = (t_float *)t_resizebytes(x->frontImag, 0, sizeof(float) * (x->fftHalfSize+1));
  //
  x->backReal = (t_float *)t_getbytes(0);
  x->backReal = (t_float *)t_resizebytes(x->backReal, 0, sizeof(float) * (x->fftHalfSize+1));
  //
  x->backImag = (t_float *)t_getbytes(0);
  x->backImag = (t_float *)t_resizebytes(x->backImag, 0, sizeof(float) * (x->fftHalfSize+1));
  //
  x->prevReal = (t_float *)t_getbytes(0);
  x->prevReal = (t_float *)t_resizebytes(x->prevReal, 0, sizeof(float) * (x->fftHalfSize+1));
  //
  x->prevImag = (t_float *)t_getbytes(0);
  x->prevImag = (t_float *)t_resizebytes(x->prevImag, 0, sizeof(float) * (x->fftHalfSize+1));
  //
  x->conjReal = (t_float *)t_getbytes(0);
  x->conjReal = (t_float *)t_resizebytes(x->conjReal, 0, sizeof(float) * (x->fftHalfSize+1));
  //
  x->conjImag = (t_float *)t_getbytes(0);
  x->conjImag = (t_float *)t_resizebytes(x->conjImag, 0, sizeof(float) * (x->fftHalfSize+1));
  //
  x->outReal = (t_float *)t_getbytes(0);
  x->outReal = (t_float *)t_resizebytes(x->outReal, 0, sizeof(float) * (x->fftHalfSize+1));
  //
  x->outImag = (t_float *)t_getbytes(0);
  x->outImag = (t_float *)t_resizebytes(x->outImag, 0, sizeof(float) * (x->fftHalfSize+1));
  //
  x->rsqrt = (t_float *)t_getbytes(0);
  x->rsqrt = (t_float *)t_resizebytes(x->rsqrt, 0, sizeof(float) * (x->fftHalfSize+1));
  //
  x->outSpectra = (t_float *)t_getbytes(0);
  x->outSpectra = (t_float *)t_resizebytes(x->outSpectra, 0, sizeof(float) * x->fftSize);
  //
  x->analWindow = (t_float *)t_getbytes(0);
  x->analWindow = (t_float *)t_resizebytes(x->analWindow, 0, sizeof(float) * x->fftSize);
  //
  x->synthWindow = (t_float *)t_getbytes(0);
  x->synthWindow = (t_float *)t_resizebytes(x->synthWindow, 0, sizeof(float) * x->fftSize);
 //
  x->nonoverlappedOutBuff = (t_float *)t_getbytes(0);
  x->nonoverlappedOutBuff = (t_float *)t_resizebytes(x->nonoverlappedOutBuff, 0, sizeof(float) * x->fftSize * x->overlap);
 

  init_window(x);



 
 

return (void *)x;
}

	/***********init window*********/
static void init_window(t_phasevoc_tilde *x)
{
  int i;
  for ( i = 0; i < x->fftSize; i++ )
    {
    //hann
     x->analWindow[i] = x->synthWindow[i] = 0.5f * (1 - cosf(x->twoPi*i/(x->fftSize-1)));
     x->synthWindow[i] *= 2.0/3.0/x->fftSize;
    }
     //hamm
  //x->analWindow[i] = x->synthWindow[i] = (float) (.54f - (.46f * cosf(x->twoPi * i / (x->fftSize - 1)) ) );
    
 
}

/****************************destructor*********************************/
static void phasevoc_tilde_free(t_phasevoc_tilde *x)
{
  if(x->circBuff != 0) t_freebytes(x->circBuff, sizeof(float) * x->circBuffLength);
  if(x->sampHoldBuff != 0) t_freebytes(x->sampHoldBuff, sizeof(float) * x->sampHoldBuffLength);
  if(x->inFrontBuff != 0) t_freebytes(x->inFrontBuff, sizeof(float) * x->fftSize);
  if(x->inBackBuff != 0) t_freebytes(x->inBackBuff, sizeof(float) * x->fftSize);
  if(x->outBuff != 0) t_freebytes(x->outBuff, sizeof(float) * x->fftSize);
  if(x->inFrontShift != 0) t_freebytes(x->inFrontShift, sizeof(float) * x->fftSize);
  if(x->inBackShift != 0) t_freebytes(x->inBackShift, sizeof(float) * x->fftSize);
  if(x->frontReal != 0) t_freebytes(x->frontReal, sizeof(float) * (x->fftHalfSize+1));
  if(x->frontImag != 0) t_freebytes(x->frontImag, sizeof(float) * (x->fftHalfSize+1));
  if(x->backReal != 0) t_freebytes(x->backReal, sizeof(float) * (x->fftHalfSize+1));
  if(x->backImag != 0) t_freebytes(x->backImag, sizeof(float) * (x->fftHalfSize+1));
  if(x->prevReal != 0) t_freebytes(x->prevReal, sizeof(float) * (x->fftHalfSize+1));
  if(x->prevImag != 0) t_freebytes(x->prevImag, sizeof(float) * (x->fftHalfSize+1));
  if(x->conjReal != 0) t_freebytes(x->conjReal, sizeof(float) * (x->fftHalfSize+1));
  if(x->conjImag != 0) t_freebytes(x->conjImag, sizeof(float) * (x->fftHalfSize+1));
  if(x->outReal != 0) t_freebytes(x->outReal, sizeof(float) * (x->fftHalfSize+1));
  if(x->outImag != 0) t_freebytes(x->outImag, sizeof(float) * (x->fftHalfSize+1));
  if(x->rsqrt != 0) t_freebytes(x->rsqrt, sizeof(float) * (x->fftHalfSize+1));
  if(x->outSpectra != 0) t_freebytes(x->outSpectra, sizeof(float) * x->fftSize);
  if(x->analWindow != 0) t_freebytes(x->analWindow, sizeof(float) * x->fftSize);
  if(x->synthWindow != 0) t_freebytes(x->synthWindow, sizeof(float) * x->fftSize);
  if(x->nonoverlappedOutBuff != 0) t_freebytes(x->nonoverlappedOutBuff, sizeof(float) * x->fftSize * x->overlap);
  

}
/**************dsp****************/
void phasevoc_tilde_dsp(t_phasevoc_tilde *x, t_signal **sp)
{
 
  dsp_add(phasevoc_tilde_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
 
}
/******************setup**********/
void phasevoc_tilde_setup(void)
{
  phasevoc_tilde_class = class_new(gensym("phasevoc~"),
				   (t_newmethod)phasevoc_tilde_new,
				   (t_method)phasevoc_tilde_free,
				   sizeof(t_phasevoc_tilde),
				   0,
				   A_GIMME,
				   0);
  
CLASS_MAINSIGNALIN(phasevoc_tilde_class, t_phasevoc_tilde, f);
  
  class_addmethod(phasevoc_tilde_class, 
		  (t_method)phasevoc_tilde_timestretch, 
		  gensym("timestretch"), 
		  A_DEFFLOAT, 
		  0);
  
  class_addmethod(phasevoc_tilde_class,
		  (t_method)phasevoc_tilde_pitchshift,
		  gensym("pitchshift"), 
		  A_DEFFLOAT, 
		  0);
  
  class_addmethod(phasevoc_tilde_class, 
		  (t_method)phasevoc_tilde_samphold, 
		  gensym("samphold"), 
		  A_DEFFLOAT, 
		  0);
  
  class_addmethod(phasevoc_tilde_class, 
		  (t_method)phasevoc_tilde_phaselock, 
		  gensym("phaselock"), 
		  A_DEFFLOAT, 
		  0);
  
  
  class_addmethod(phasevoc_tilde_class, 
		  (t_method)phasevoc_tilde_reinit,
		  gensym("reinit"), 
		  0);


  class_addmethod(phasevoc_tilde_class, 
		  (t_method)phasevoc_tilde_pause, 
		  gensym("pause"),
		  0);  

  class_addmethod(phasevoc_tilde_class, 
		  (t_method)phasevoc_tilde_unpause,
		  gensym("unpause"), 
		  0);

  
  class_addmethod(phasevoc_tilde_class, 
		  (t_method)phasevoc_tilde_dsp, 
		  gensym("dsp"), 
		  0);

  //  class_sethelpsymbol(phasevoc_tilde_class, gensym("help-phasevoc~"));

}


