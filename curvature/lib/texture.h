/***************************************************************************
texture.h
----------------------------------------------------------------------------
begin                : june 2003
copyright            : (C) 2003 by Pierre Alliez - INRIA
email                : pierre.alliez@sophia.inria.fr
***************************************************************************/

#ifndef _TEXTURE_
#define _TEXTURE_


class CTexture
{

// Member data
private :

  unsigned char *m_pData;         // data
  unsigned int   m_Width;         // width (pixels)
  unsigned int   m_Height;        // height (pixels)
  unsigned int   m_Depth;         // bits per pixel 
  char          *m_pFileName;     // texture image file name

  BITMAPINFOHEADER m_Header;      // image header (display on device context)
  unsigned int     m_WidthByte32; // width (in bytes, and 32 bits aligned)

public :

  // Construction / destruction
  CTexture();
  virtual ~CTexture();

  // File reading
  int ReadFile(char *filename,unsigned int width=-1,
	       unsigned int height=-1,unsigned int depth=-1);
  int ReadFileBMP(char *filename);
  int ReadFileRAW(char *filename,unsigned int width,
		  unsigned int height,unsigned int depth);

  // File saving
  int SaveFile(char *filename);
  int SaveFileBMP(char *filename);
  int SaveFileRAW(char *filename);
	
  // Get Data (explicit inline functions)
  unsigned char *GetData()     { return m_pData; }
  unsigned int   GetWidth()    { return m_Width; }
  unsigned int   GetWidthByte32()    { return m_WidthByte32; }
  unsigned int   GetHeight()   { return m_Height;}
  unsigned int   GetDepth()    { return m_Depth; }
  char*          GetFileName() { return m_pFileName; }

  // Misc
  int IsValid();
  int SameSize(CTexture *pTexture);
  int BGRtoRGB();
  static int HigherPowerOfTwo(int value);
  static int LowerPowerOfTwo(int value);

  // Updating
  void UpdateWidthByte32();
  void UpdateHeader();
  unsigned int WidthByte32(unsigned int width,unsigned int depth);

  // Clipboard
  HANDLE ExportHandle();
	int ImportHandle(HANDLE handle);

  // Alpha
  int HasAlpha() { return (m_Depth == 32); }
  int AddAlphaLayer(unsigned char alpha);
  int SetAlphaLayer(unsigned char alpha);
  int SetAlphaLayer(CTexture *pTexture);	// Put an alpha layer from grey scale

  // DuplicateMirror
  int DuplicateMirror(int left=0,int top=0,int right=-1,int bottom=-1);
  int DuplicateRepeatWidth(int left=0,int top=0,int right=-1,int bottom=-1);
  int Extract(int left=0,int top=0,int right=-1,int bottom=-1);

  // Display
  int Draw(CDC *pDC,int xOffset=0,int yOffset=0, int width=-1, int height=-1);
	int Stretch(CDC *pDC,CRect *pRect);

  // Buffer
  int ReadBuffer(unsigned char *buffer, int width, int height, int depth);
	int ReadBufferByte32(unsigned char *pData,int width,int height);
  int ReadBuffer(float *buffer, int width, int height, int depth);
	int	ReadBuffer(double *buffer,int width, int height, int depth);
  int  Grey(unsigned int x,unsigned int y);
  void Color(unsigned int x,unsigned int y,
	     unsigned char *pRed,unsigned char *pGreen,unsigned char *pBlue);
	int ReadBuffer(float **ppBuffer, int width, int height,float ratio = 1.0f);
	int WriteBuffer(float **ppBuffer,int width,int height);
	int WriteBuffer32(float **ppBuffer,int width,int height);

	// Misc
	void GenerateGrid(unsigned int width,unsigned int height,int size,
		unsigned char r,unsigned char g,unsigned char b,
		unsigned char rb,unsigned char gb,unsigned char bb);
	void GenerateLines(unsigned int width,unsigned int height,int size,
		unsigned char r,unsigned char g,unsigned char b,
		unsigned char rb,unsigned char gb,unsigned char bb);

  // Memory
  int  Alloc(unsigned int width,unsigned int height,unsigned int depth);
  void Free();
	void Fill(unsigned char r,unsigned char g,unsigned char b);
	void Copy(CTexture *pTexture);

	void GreyToColor(unsigned char grey,unsigned char r, 
		unsigned char g,unsigned char b);
	void ColorToColor(unsigned char r1,unsigned char g1,
		unsigned char b1,unsigned char r2, unsigned char g2,
		unsigned char b2);
};

#endif // _TEXTURE_
