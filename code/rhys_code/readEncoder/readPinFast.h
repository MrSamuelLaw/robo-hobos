bool readPinFast(uint8_t pin)
{
  switch(pin)
  {
//    case 0:
//      return (*( (volatile uint8_t*) 0x2C ) & (uint8_t)1) != 0;
//    case 1:
//      return (*( (volatile uint8_t*) 0x2C ) & (uint8_t)2) != 0;
    case 2:
      return (*( (volatile uint8_t*) 0x2C ) & (uint8_t)16) != 0;
//    case 3:
//      return (*( (volatile uint8_t*) 0x2C ) & (uint8_t)32) != 0;
    case 4:
      return (*( (volatile uint8_t*) 0x32 ) & (uint8_t)32) != 0;
    case 5:
      return (*( (volatile uint8_t*) 0x2C ) & (uint8_t)8) != 0;
    case 6:
      return (*( (volatile uint8_t*) 0x100 ) & (uint8_t)8) != 0;
    case 7:
      return (*( (volatile uint8_t*) 0x100 ) & (uint8_t)16) != 0;
//    case 8:
//      return (*( (volatile uint8_t*) 0x100 ) & (uint8_t)32) != 0;
//    case 9:
//      return (*( (volatile uint8_t*) 0x100 ) & (uint8_t)64) != 0;
    case 10:
      return (*( (volatile uint8_t*) 0x23 ) & (uint8_t)16) != 0;
//    case 11:
//      return (*( (volatile uint8_t*) 0x23 ) & (uint8_t)32) != 0;
//    case 12:
//      return (*( (volatile uint8_t*) 0x23 ) & (uint8_t)64) != 0;
//    case 13:
//      return (*( (volatile uint8_t*) 0x23 ) & (uint8_t)128) != 0;
    case 14:
      return (*( (volatile uint8_t*) 0x103 ) & (uint8_t)2) != 0;
    case 15:
      return (*( (volatile uint8_t*) 0x103 ) & (uint8_t)1) != 0;
    case 16:
      return (*( (volatile uint8_t*) 0x100 ) & (uint8_t)2) != 0;
    case 17:
      return (*( (volatile uint8_t*) 0x100 ) & (uint8_t)1) != 0;
    case 18:
      return (*( (volatile uint8_t*) 0x29 ) & (uint8_t)8) != 0;
    case 19:
      return (*( (volatile uint8_t*) 0x29 ) & (uint8_t)4) != 0;
    case 20:
      return (*( (volatile uint8_t*) 0x29 ) & (uint8_t)2) != 0;
    case 21:
      return (*( (volatile uint8_t*) 0x29 ) & (uint8_t)1) != 0;
    case 22:
      return (*( (volatile uint8_t*) 0x20 ) & (uint8_t)1) != 0;
    case 23:
      return (*( (volatile uint8_t*) 0x20 ) & (uint8_t)2) != 0;
    case 24:
      return (*( (volatile uint8_t*) 0x20 ) & (uint8_t)4) != 0;
    case 25:
      return (*( (volatile uint8_t*) 0x20 ) & (uint8_t)8) != 0;
    case 26:
      return (*( (volatile uint8_t*) 0x20 ) & (uint8_t)16) != 0;
    case 27:
      return (*( (volatile uint8_t*) 0x20 ) & (uint8_t)32) != 0;
    case 28:
      return (*( (volatile uint8_t*) 0x20 ) & (uint8_t)64) != 0;
    case 29:
      return (*( (volatile uint8_t*) 0x20 ) & (uint8_t)128) != 0;
    case 30:
      return (*( (volatile uint8_t*) 0x26 ) & (uint8_t)128) != 0;
    case 31:
      return (*( (volatile uint8_t*) 0x26 ) & (uint8_t)64) != 0;
    case 32:
      return (*( (volatile uint8_t*) 0x26 ) & (uint8_t)32) != 0;
    case 33:
      return (*( (volatile uint8_t*) 0x26 ) & (uint8_t)16) != 0;
    case 34:
      return (*( (volatile uint8_t*) 0x26 ) & (uint8_t)8) != 0;
    case 35:
      return (*( (volatile uint8_t*) 0x26 ) & (uint8_t)4) != 0;
    case 36:
      return (*( (volatile uint8_t*) 0x26 ) & (uint8_t)2) != 0;
    case 37:
      return (*( (volatile uint8_t*) 0x26 ) & (uint8_t)1) != 0;
    case 38:
      return (*( (volatile uint8_t*) 0x29 ) & (uint8_t)128) != 0;
    case 39:
      return (*( (volatile uint8_t*) 0x32 ) & (uint8_t)4) != 0;
    case 40:
      return (*( (volatile uint8_t*) 0x32 ) & (uint8_t)2) != 0;
    case 41:
      return (*( (volatile uint8_t*) 0x32 ) & (uint8_t)1) != 0;
    case 42:
      return (*( (volatile uint8_t*) 0x109 ) & (uint8_t)128) != 0;
    case 43:
      return (*( (volatile uint8_t*) 0x109 ) & (uint8_t)64) != 0;
    case 44:
      return (*( (volatile uint8_t*) 0x109 ) & (uint8_t)32) != 0;
    case 45:
      return (*( (volatile uint8_t*) 0x109 ) & (uint8_t)16) != 0;
    case 46:
      return (*( (volatile uint8_t*) 0x109 ) & (uint8_t)8) != 0;
    case 47:
      return (*( (volatile uint8_t*) 0x109 ) & (uint8_t)4) != 0;
    case 48:
      return (*( (volatile uint8_t*) 0x109 ) & (uint8_t)2) != 0;
    case 49:
      return (*( (volatile uint8_t*) 0x109 ) & (uint8_t)1) != 0;
    case 50:
      return (*( (volatile uint8_t*) 0x23 ) & (uint8_t)8) != 0;
    case 51:
      return (*( (volatile uint8_t*) 0x23 ) & (uint8_t)4) != 0;
    case 52:
      return (*( (volatile uint8_t*) 0x23 ) & (uint8_t)2) != 0;
    case 53:
      return (*( (volatile uint8_t*) 0x23 ) & (uint8_t)1) != 0;
    default:
      return 0;
  }
}
