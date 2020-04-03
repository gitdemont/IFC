#include <Rcpp.h>
#include <progress.hpp>
#include <iostream>
#include <fstream>
#include "../inst/include/utils.hpp"
#include "../inst/include/color.hpp"
#include "../inst/include/trans.hpp"
#include "../inst/include/matrix_logic.hpp"
#include "../inst/include/resize.hpp"
using namespace Rcpp;

// template to swap bytes for each types
template <typename T> T bytes_swap(T val) {
  T out;
  char *p_val = (char*)&val;
  char *p_out = (char*)&out;
  int size = sizeof(T);
  for(int i=0; i<size; i++) {
    p_out[size-1-i] = p_val[i];
  }
  return out;
}

//' @title File Scanner
//' @name cpp_scanFirst
//' @description
//' Scans file for 1st occurence of a target string.
//' If found, it returns the position in bytes of the target.
//' Otherwise, it returns 0.
//' @param fname string, path to file.
//' @param target string, exact string to be searched for. At least 1 character and should not exceed 1024 characters.
//' @param start size_t, position where to begin search.
//' It can't be superior or equal than file size or end (when end is different from 0 and inferior than file size).
//' @param end size_t, position where to stop searching. Default is 0.
//' Search will end up at this position unless it is higher than file size.
//' In such case, search will end up when file end will be reached.
//' @param buf_size uint8_t, size of buffer used to search for target (in kilo-Bytes, will be forced to be between 2 and 1024). Default is 64.
//' @return size_t index of first target character found within target plus 1 or 0 if not found.
//' @keywords internal
////' @export
// [[Rcpp::export]]
size_t cpp_scanFirst(const std::string fname, const std::string target, const size_t start = 0, const size_t end = 0, const uint8_t buf_size = 64) {
  uint16_t L = target.length();
  if(L < 1) {
    Rcerr <<  "cpp_scanFirst: target should be at least 1 character";
    Rcpp::stop("cpp_scanFirst: target should be at least 1 character");
  }
  if(L > 1024) {
    Rcerr <<  "cpp_scanFirst: target should not exceed 1024 characters but is of length: [" << L << "], for file:" << std::endl << fname << std::endl;
    Rcpp::stop("cpp_scanFirst: target should not exceed 1024 characters");
  }
  
  std::ifstream fi(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
  if (fi.is_open()) {
    try {
      fi.seekg(0, std::ios::end);
      size_t filesize = fi.tellg(), end_at = 0, n = 0, pos = 0;
      bool keep_searching = true;
      
      // determines where to stop
      end_at = filesize;
      if(end != 0) if(end < filesize) end_at = end;
      if(start > (end_at - L)) return 0;

      // clip buffer to [2,255]
      uint8_t buf_s = buf_size;
      if(buf_s < 2) buf_s = 2;
      uint32_t s = 1024 * buf_s; // allow reading chunk of size buf_s * kilo-Bytes
      
      fi.seekg(start, std::ios::beg);
      while(keep_searching) {
        pos = fi.tellg();
        if(pos >= end_at) break; // break loop when end of the file is reached (or end when provided);
        if(pos > (start + L)) pos -= L; // make a frame shift for every chunk except 1st one to be sure to include target when it is positioned between 2 chunks;
        fi.seekg(pos, std::ios::beg);
        if((end_at - pos) < s) s = end_at - pos; // reduce buffer size for last chunk to avoid reading outside of the file;
        std::vector<char> buffer(s);
        fi.read(buffer.data(), s);
        std::string str(buffer.begin(), buffer.end());
        n = str.find(target);
        if(n != std::string::npos) keep_searching = false;
      }
      fi.close();
      if(!keep_searching) return pos + n + 1;
    }
    catch(std::exception &ex) { 
      fi.close();
      forward_exception_to_r(ex);
    }
    catch(...) {
      Rcpp::stop("cpp_scanFirst: c++ exception (unknown reason)");
    }
  }
  else {
    Rcerr << "cpp_scanFirst: Unable to open " << fname << std::endl;
    Rcpp::stop("cpp_scanFirst: Unable to open file");
  }
  return 0;
}

//' @title TIFF Checker
//' @name cpp_checkTIFF
//' @description
//' Checks if file is a TIFF.
//' @details If file is a TIFF it returns endianness of file, 'big' or 'little.
//' Otherwise, it shows an error and returns an empty string.
//' @param fname string, path to file.
//' @source TIFF 6.0 specifications available at \url{https://www.adobe.io/open/standards/TIFF.html}
//' @keywords internal
////' @export
// [[Rcpp::export]]
std::string cpp_checkTIFF (std::string fname) {
  std::ifstream fi(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
  std::string out = "";
  if (fi.is_open()) {
    try {
      fi.seekg(0, std::ios::end);
      size_t filesize = fi.tellg();
      if(filesize < 22) { // the 4 magic bytes + 1st Offset (4 bytes) + Count of Entry (2 bytes) + at least 1 IFD Entry (12 bytes) + last offset (4 bytes)
        Rcpp::stop("cpp_checkTIFF: File is too small");
      }
      uint16_t magic;
      fi.seekg(0, std::ios::beg);
      fi.read((char *)&magic,sizeof(magic));
      
      if(magic == 18761) out = "little"; // 49,49
      if(magic == 19789) out = "big"; // 4D,4D
      if(out == "") {
        Rcerr <<  "cpp_checkTIFF: " << fname << "\nis not a XIF file: No magic bytes 0-1" << std::endl;
        Rcpp::stop("cpp_checkTIFF: File is not a XIF file: No magic bytes 0-1");
      }
      fi.read((char *)&magic,sizeof(magic));
      if(out == "big") magic = bytes_swap(magic);
      if(magic != 42) {
        Rcerr <<  "cpp_checkTIFF: " << fname << "\nis not a XIF file: No magic bytes 2-3" << std::endl;
        Rcpp::stop("cpp_checkTIFF: File is not a XIF file: No magic bytes 2-3");
      }
      fi.close();
      return out;
    }
    catch(std::exception &ex) {	
      fi.close();
      forward_exception_to_r(ex);
    }
    catch(...) { 
      Rcpp::stop("cpp_checkTIFF: c++ exception (unknown reason)"); 
    }
  }
  else {
    Rcerr << "cpp_checkTIFF: Unable to open " << fname << std::endl;
    Rcpp::stop("cpp_checkTIFF: Unable to open file");
  }
  return "";
}

//' @title IFC_offsets Computation without Id Determination
//' @name cpp_getoffsets_noid
//' @description
//' Returns offsets of the IFDs (Image Field Directory) within a TIFF file.
//' @param fname string, path to file.
//' @param obj_count R_len_t, numbers of objects present in the file. Default is 0.
//' If obj_count <= 0 then progress_bar is forced to false.
//' @param display_progress bool, whether to display a progress bar. Default is false.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @source TIFF 6.0 specifications available at \url{https://www.adobe.io/open/standards/TIFF.html}
//' @return an integer vector with offsets of IFDs found.
//' @keywords internal
////' @export
// [[Rcpp::export]]
IntegerVector cpp_getoffsets_noid(std::string fname, R_len_t obj_count = 0, bool display_progress = false, bool verbose = false) {
  bool swap = false;
  std::string endianness = cpp_checkTIFF(fname);
  if(endianness == "big") swap = true;
  
  std::ifstream fi(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
  if (fi.is_open()) {
    try {
      fi.seekg(0, std::ios::end);
      size_t filesize = fi.tellg(); // offsets are uint32 so we can't go further
      uint32_t offset = 4; // offsets are uint32 
      size_t pos;
      IntegerVector out;
      if(obj_count <= 0) display_progress = false;
      char buf_entries [2];
      char buf_offset [4];
      uint16_t entries; // entries are uint16

      if(verbose) Rcout << "Extracting offsets from " << fname << std::endl;
      fi.seekg(offset, std::ios::beg);
      fi.read((char*)&buf_offset, sizeof(buf_offset));
      memcpy(&offset, buf_offset, sizeof(offset));

      Progress p(obj_count * 2 + 1, display_progress);
      if(swap) {
        offset = bytes_swap(offset);
        if(!offset) {
          Rcerr << "cpp_getoffsets_noid: No IFD offsets found in\n" << fname << std::endl;
          Rcpp::stop("cpp_getoffsets_noid: No IFD offsets found");
        }
        while(offset) {
          p.increment();
          out.push_back(offset);
          fi.seekg(offset, std::ios::beg);
          fi.read((char*)buf_entries, sizeof(buf_entries));
          memcpy(&entries, buf_entries, sizeof(entries));
          entries = bytes_swap(entries);
          pos = 2 + offset + 12 * entries;
          if((pos > filesize) || (pos >= 0xffffffff)) { // can't read to more than 4,294,967,295
            Rcerr << "cpp_getoffsets_noid: Buffer overrun in\n" << fname << std::endl;
            Rcpp::stop("cpp_getoffsets_noid: Buffer overrun");
          }
          fi.seekg(pos, std::ios::beg);
          fi.read((char*)buf_offset, sizeof(buf_offset));
          memcpy(&offset, buf_offset, sizeof(offset));
          offset = bytes_swap(offset);
          if(Progress::check_abort()) {
            Rcerr << "cpp_getoffsets_noid: Interrupted by user" << std::endl;
            Rcpp::stop("cpp_getoffsets_noid: Interrupted by user");
          }
        }
      } else {
        if(!offset) {
          Rcerr << "cpp_getoffsets_noid: No IFD offsets found in\n" << fname << std::endl;
          Rcpp::stop("cpp_getoffsets_noid: No IFD offsets found");
        }
        while(offset) {
          p.increment();
          out.push_back(offset);
          fi.seekg(offset, std::ios::beg);
          fi.read((char*)buf_entries, sizeof(buf_entries));
          memcpy(&entries, buf_entries, sizeof(entries));
          pos = 2 + offset + 12 * entries;
          if((pos > filesize) || (pos >= 0xffffffff)) { // can't read to more than 4,294,967,295
            Rcerr << "cpp_getoffsets_noid: Buffer overrun in\n" << fname << std::endl;
            Rcpp::stop("cpp_getoffsets_noid: Buffer overrun");
          }
          fi.seekg(pos, std::ios::beg);
          fi.read((char*)buf_offset, sizeof(buf_offset));
          memcpy(&offset, buf_offset, sizeof(offset));
          if(Progress::check_abort()) {
            Rcerr << "cpp_getoffsets_noid: Interrupted by user" << std::endl;
            Rcpp::stop("cpp_getoffsets_noid: Interrupted by user");
          }
        }
      }
      fi.close();
      return out;
    }
    catch(std::exception &ex) {	
      fi.close();
      forward_exception_to_r(ex);
    }
    catch(...) { 
      Rcpp::stop("cpp_getoffsets_noid: c++ exception (unknown reason)"); 
    }
  }
  else {
    Rcerr << "cpp_getoffsets_noid: Unable to open " << fname << std::endl;
    Rcpp::stop("cpp_getoffsets_noid: Unable to open file");
  }
  return IntegerVector::create(0);
}

//' @title IFD Tags Extraction
//' @name cpp_getTAGS
//' @description
//' Returns TAGS contained within an IFD (Image Field Directory) entry.
//' @param fname string, path to file.
//' @param offset uint32_t, position of the IFD beginning.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is 'false'.
//' @param trunc_bytes uint8_t, number of bytes to extract for bytes/strings TAGS (1, 2, 6 or 7). Default is 22.
//' @param force_trunc bool, whether to force truncation for all TAGS types. Default is 'false'.
//' @source TIFF 6.0 specifications available at \url{https://www.adobe.io/open/standards/TIFF.html}
//' @keywords internal
////' @export
// [[Rcpp::export]]
List cpp_getTAGS (std::string fname, uint32_t offset, bool verbose = false, uint8_t trunc_bytes = 22, bool force_trunc = false) {
  bool swap = false;
  std::string endianness = cpp_checkTIFF(fname); // it may be better to remove this check each time
  if(endianness == "big") swap = true; 
  
  std::ifstream fi(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
  if (fi.is_open()) {
    try {
      fi.seekg(0, std::ios::end);
      size_t filesize = fi.tellg();
      fi.seekg(0, std::ios::beg);
      if(offset > (filesize - 2)) {
        Rcerr << "cpp_getTAGS: offset@" << offset << " points to outside of " << fname << std::endl;
        Rcpp::stop("cpp_getTAGS: TAG offset is higher than file size");
      }
      if(verbose) Rcout << "Extracting TAGs from "<< fname << ", filesize:" << filesize << " @offset:" << offset << std::endl;
      
      uint16_t entries;
      uint16_t IFD_tag;
      uint16_t IFD_type;
      uint32_t IFD_count;
      uint32_t IFD_value;
      uint32_t IFD_bytes;
      uint32_t i;
      uint32_t L;
      bool IFD_off;
      bool is_char;
      size_t pos;
      // R_len_t i_map;
      char buf_entries [2];
      char buf_dir_entry [12];
      
      IntegerVector sizes = IntegerVector::create(0,1,1,2,4,4,1,1,2,4,4,4,8);
      IntegerVector multi = IntegerVector::create(0,1,1,1,1,2,1,1,1,1,2,1,1);
      
      fi.seekg(offset, std::ios::beg);
      fi.read((char*)&buf_entries, sizeof(buf_entries));
      memcpy(&entries, buf_entries, sizeof(entries));
      if(swap) entries = bytes_swap(entries);
      pos = offset + sizeof(buf_entries);
      
      RObject IMAGE_LENGTH = R_NilValue;          // 257
      RObject IMAGE_WIDTH = R_NilValue;           // 256
      // RObject TILE_LENGTH = R_NilValue;        // 323
      // RObject TILE_WIDTH = R_NilValue;         // 322
      RObject OBJECT_ID = R_NilValue;             // 33003
      // IntegerVector OBJECT_ID_2 = IR_NilValue; // 33024
      RObject COMPRESSION = R_NilValue;           // 259
      RObject TYPE = R_NilValue;                  // 33002
      RObject STRIP_OFFSETS = R_NilValue;         // 273
      RObject STRIP_BYTE_COUNTS = R_NilValue;     // 279
      // RObject TILE_OFFSETS = R_NilValue;       // 324
      // RObject TILE_BYTE_COUNTS = IR_NilValue;  // 325
      RObject BG_MEAN = R_NilValue;               // 33052
      RObject BG_STD = R_NilValue;                // 33053

      List INFOS(entries);
      CharacterVector NAMES(entries);
      if(verbose) Rcout << "Entries: " << entries << std::endl;
      if(swap) {
        for(uint16_t k = 0; k < entries; k++) {
          if((pos > (filesize - 12)) || (pos >= 0xffffffff)) {
            Rcerr << "cpp_getTAGS: in IFD: " << k << " @" << pos + 12 << " is outside of " << fname << std::endl;
            Rcpp::stop("cpp_getTAGS: IFD detected position is higher than file size");
          }
          fi.seekg(pos, std::ios::beg);
          fi.read((char*)&buf_dir_entry, sizeof(buf_dir_entry));
          memcpy(&IFD_tag, buf_dir_entry, sizeof(IFD_tag));
          IFD_tag = bytes_swap(IFD_tag);
          memcpy(&IFD_type, buf_dir_entry + 2, sizeof(IFD_type));
          IFD_type = bytes_swap(IFD_type);
          if((IFD_type > 12) || (IFD_type < 1)) {
            Rcerr << "cpp_getTAGS: " << fname << std::endl;
            Rcerr << "cpp_getTAGS: in IFD: " << k << " IFD_type=" << IFD_type << " is not allowed" << std::endl;
            Rcpp::stop("cpp_getTAGS: Value not allowed for IFD type");
          }
          memcpy(&IFD_count, buf_dir_entry + 4, sizeof(IFD_count));
          IFD_count = bytes_swap(IFD_count);
          memcpy(&IFD_value, buf_dir_entry + 8, sizeof(IFD_value));
          IFD_value = bytes_swap(IFD_value);
          pos = fi.tellg();
          
          if((IFD_count == 0) && (IFD_value != 0)) IFD_count = 1; // fix for compatibility
          IFD_bytes = sizes[IFD_type] * multi[IFD_type] * IFD_count;
          L = multi[IFD_type] * IFD_count;
          is_char = (IFD_type == 1) || (IFD_type == 2) || (IFD_type == 6) || (IFD_type == 7);
          if(is_char || force_trunc) {
            if(L > trunc_bytes) {
              L = trunc_bytes;
            }
          }
          
          if(IFD_bytes > 4) {
            IFD_off = true;
            if((IFD_value + IFD_bytes) > filesize) {
              Rcerr << "cpp_getTAGS: in IFD: " << k << " @" << IFD_value + IFD_bytes << " is outside of " << fname << std::endl;
              Rcpp::stop("cpp_getTAGS: IFD value points to outside of file");
            }
            fi.seekg(IFD_value, std::ios::beg);
          } else {
            fi.seekg(pos - 4, std::ios::beg);
            IFD_off = false;
          }
          if(verbose) Rcout << "Tag:" << IFD_tag << " Typ:" << IFD_type << " Count:" << IFD_count << " Value:" << IFD_value << " Bytes:" << IFD_bytes << " Off:" << IFD_off << std::endl;
          NAMES[k] = to_string(IFD_tag);
          
          if(is_char) {
            if(IFD_type == 2) {
              std::vector<char> buf_offset(L * sizes[IFD_type]);
              fi.read(buf_offset.data(), L * sizes[IFD_type]);
              std::string IFD_map(buf_offset.begin(), buf_offset.end());
              INFOS[k] = List::create(_["tag"] = IFD_tag,
                                      _["typ"] = IFD_type, 
                                      _["siz"] = IFD_count, 
                                      _["val"] = IFD_value, 
                                      _["byt"] = IFD_bytes, 
                                      _["len"] = IFD_count * multi[IFD_type], 
                                      _["off"] = IFD_off,
                                      _["map"] = IFD_map);
            } else {
              RawVector IFD_map(L);
              switch(IFD_type) {
              case 1: {
                unsigned char ele;
                for(i = 0; i < L; i++) {
                  fi.read((char *)&ele, sizes[IFD_type]);
                  IFD_map[i] = ele;
                }
              }
                break;
              case 6: {
                signed char ele;
                for(i = 0; i < L; i++) {
                  fi.read((char *)&ele, sizes[IFD_type]);
                  IFD_map[i] = ele;
                }
              }
                break;
              case 7: {
                char ele;
                for(i = 0; i < L; i++) {
                  fi.read((char *)&ele, sizes[IFD_type]);
                  IFD_map[i] = ele;
                }
              }
                break;
              }
              INFOS[k] = List::create(_["tag"] = IFD_tag,
                                      _["typ"] = IFD_type, 
                                      _["siz"] = IFD_count, 
                                      _["val"] = IFD_value, 
                                      _["byt"] = IFD_bytes, 
                                      _["len"] = IFD_count * multi[IFD_type], 
                                      _["off"] = IFD_off,
                                      _["map"] = IFD_map);
            }
          } else {
            NumericVector IFD_map(L);
            switch(IFD_type) {
            case 3: {
              uint16_t ele;
              for(i = 0; i < L; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = bytes_swap(ele);
              }
            }
              break;
            case 4: {
              uint32_t ele;
              for(i = 0; i < L; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = bytes_swap(ele);
              }
            }
              break;
            case 5: {
              uint32_t ele;
              for(i = 0; i < L; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = bytes_swap(ele);
              }
            }
              break;
            case 8: {
              int16_t ele;
              for(i = 0; i < L; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = bytes_swap(ele);
              }
            }
              break;
            case 9: {
              int32_t ele;
              for(i = 0; i < L; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = bytes_swap(ele);
              }
            }
              break;
            case 10: {
              int32_t ele;
              for(i = 0; i < L; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = bytes_swap(ele);
              }
            }
              break;
            case 11: {
              float ele;
              for(i = 0; i < L; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = bytes_swap(ele);
              }
            }
              break;
            case 12: {
              double ele;
              for(i = 0; i < L; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = bytes_swap(ele);
              }
            }
              break;
            }
            
            switch(IFD_tag) {
            // TODO, ask AMNIS whether they use the TILE (OFFSETS/BYTE_COUNTS) or STRIP (OFFSETS/BYTE_COUNTS)
            // TODO, ask AMNIS whether they use the TILE (LENGTH/WIDTH) or IMAGE (LENGTH/WIDTH)
            case 256: {
              IMAGE_WIDTH = IFD_map[0];
              break;
            }
            case 257: {
              IMAGE_LENGTH = IFD_map[0];
              break;
            }
            case 259: {
              COMPRESSION = IFD_map[0];
              break;
            }
            case 273: {
              STRIP_OFFSETS = IFD_map[0];
              break;
            }
            case 279: {
              STRIP_BYTE_COUNTS = IFD_map[0];
              break;
            }
            case 33002: {
              TYPE = IFD_map[0];
              break;
            }
            case 33003: {
              OBJECT_ID = IFD_map[0];
              break;
            }
            case 33052: {
              BG_MEAN = IFD_map;
              break;
            }
            case 33053: {
              BG_STD = IFD_map;
              break;
            }
            }
            INFOS[k] = List::create(_["tag"] = IFD_tag,
                                    _["typ"] = IFD_type, 
                                    _["siz"] = IFD_count, 
                                    _["val"] = IFD_value, 
                                    _["byt"] = IFD_bytes, 
                                    _["len"] = IFD_count * multi[IFD_type], 
                                    _["off"] = IFD_off,
                                    _["map"] = IFD_map);
          }
        }
      } else {
        for(uint16_t k = 0; k < entries; k++) {
          if((pos > (filesize - 12)) || (pos >= 0xffffffff)) {
            Rcerr << "cpp_getTAGS: in IFD: " << k << " @" << pos + 12 << " is outside of " << fname << std::endl;
            Rcpp::stop("cpp_getTAGS: IFD detected position is higher than file size");
          }
          fi.seekg(pos, std::ios::beg);
          fi.read((char*)&buf_dir_entry, sizeof(buf_dir_entry));
          memcpy(&IFD_tag, buf_dir_entry, sizeof(IFD_tag));
          memcpy(&IFD_type, buf_dir_entry + 2, sizeof(IFD_type));
          if((IFD_type > 12) || (IFD_type < 1)) {
            Rcerr << "cpp_getTAGS: " << fname << std::endl;
            Rcerr << "cpp_getTAGS: in IFD: " << k << " IFD_type=" << IFD_type << " is not allowed" << std::endl;
            Rcpp::stop("cpp_getTAGS: Value not allowed for IFD type");
          }
          memcpy(&IFD_count, buf_dir_entry + 4, sizeof(IFD_count));
          memcpy(&IFD_value, buf_dir_entry + 8, sizeof(IFD_value));
          pos = fi.tellg();
          
          if((IFD_count == 0) && (IFD_value != 0)) IFD_count = 1; // fix for compatibility
          IFD_bytes = sizes[IFD_type] * multi[IFD_type] * IFD_count;
          L = multi[IFD_type] * IFD_count;
          is_char = (IFD_type == 1) || (IFD_type == 2) || (IFD_type == 6) || (IFD_type == 7);
          if(is_char || force_trunc) {
            if(L > trunc_bytes) {
              L = trunc_bytes;
            }
          }
          
          if(IFD_bytes > 4) {
            IFD_off = true;
            if((IFD_value + IFD_bytes) > filesize) {
              Rcerr << "cpp_getTAGS: in IFD: " << k << " @" << IFD_value + IFD_bytes << " is outside of " << fname << std::endl;
              Rcpp::stop("cpp_getTAGS: IFD value points to outside of file");
            }
            fi.seekg(IFD_value, std::ios::beg);
          } else {
            fi.seekg(pos - 4, std::ios::beg);
            IFD_off = false;
          }
          if(verbose) Rcout << "Tag:" << IFD_tag << " Typ:" << IFD_type << " Count:" << IFD_count << " Value:" << IFD_value << " Bytes:" << IFD_bytes << " Off:" << IFD_off << std::endl;
          NAMES[k] = to_string(IFD_tag);
          
          if(is_char) {
            if(IFD_type == 2) {
              std::vector<char> buf_offset(L * sizes[IFD_type]);
              fi.read(buf_offset.data(), L * sizes[IFD_type]);
              std::string IFD_map(buf_offset.begin(), buf_offset.end());
              INFOS[k] = List::create(_["tag"] = IFD_tag,
                                      _["typ"] = IFD_type, 
                                      _["siz"] = IFD_count, 
                                      _["val"] = IFD_value, 
                                      _["byt"] = IFD_bytes, 
                                      _["len"] = IFD_count * multi[IFD_type], 
                                      _["off"] = IFD_off,
                                      _["map"] = IFD_map);
            } else {
              RawVector IFD_map(L);
              switch(IFD_type) {
              case 1: {
                unsigned char ele;
                for(i = 0; i < L; i++) {
                  fi.read((char *)&ele, sizes[IFD_type]);
                  IFD_map[i] = ele;
                }
              }
                break;
              case 6: {
                signed char ele;
                for(i = 0; i < L; i++) {
                  fi.read((char *)&ele, sizes[IFD_type]);
                  IFD_map[i] = ele;
                }
              }
                break;
              case 7: {
                char ele;
                for(i = 0; i < L; i++) {
                  fi.read((char *)&ele, sizes[IFD_type]);
                  IFD_map[i] = ele;
                }
              }
                break;
              }
              INFOS[k] = List::create(_["tag"] = IFD_tag,
                                      _["typ"] = IFD_type, 
                                      _["siz"] = IFD_count, 
                                      _["val"] = IFD_value, 
                                      _["byt"] = IFD_bytes, 
                                      _["len"] = IFD_count * multi[IFD_type], 
                                      _["off"] = IFD_off,
                                      _["map"] = IFD_map);
            }
          } else {
            NumericVector IFD_map(L);
            switch(IFD_type) {
            case 3: {
              uint16_t ele;
              for(i = 0; i < L; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = ele;
              }
            }
              break;
            case 4: {
              uint32_t ele;
              for(i = 0; i < L; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = ele;
              }
            }
              break;
            case 5: {
              uint32_t ele;
              for(i = 0; i < L; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = ele;
              }
            }
              break;
            case 8: {
              int16_t ele;
              for(i = 0; i < L; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = ele;
              }
            }
              break;
            case 9: {
              int32_t ele;
              for(i = 0; i < L; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = ele;
              }
            }
              break;
            case 10: {
              int32_t ele;
              for(i = 0; i < L; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = ele;
              }
            }
              break;
            case 11: {
              float ele;
              for(i = 0; i < L; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = ele;
              }
            }
              break;
            case 12: {
              double ele;
              for(i = 0; i < L; i++) {
                fi.read((char *)&ele, sizes[IFD_type]);
                IFD_map[i] = ele;
              }
            }
              break;
            }
            
            switch(IFD_tag) {
            // TODO, ask AMNIS whether they use the TILE (OFFSETS/BYTE_COUNTS) or STRIP (OFFSETS/BYTE_COUNTS)
            // TODO, ask AMNIS whether they use the TILE (LENGTH/WIDTH) or IMAGE (LENGTH/WIDTH)
            case 256: {
              IMAGE_WIDTH = IFD_map[0];
              break;
            }
            case 257: {
              IMAGE_LENGTH = IFD_map[0];
              break;
            }
            case 259: {
              COMPRESSION = IFD_map[0];
              break;
            }
            case 273: {
              STRIP_OFFSETS = IFD_map[0];
              break;
            }
            case 279: {
              STRIP_BYTE_COUNTS = IFD_map[0];
              break;
            }
            case 33002: {
              TYPE = IFD_map[0];
              break;
            }
            case 33003: {
              OBJECT_ID = IFD_map[0];
              break;
            }
            case 33052: {
              BG_MEAN = IFD_map;
              break;
            }
            case 33053: {
              BG_STD = IFD_map;
              break;
            }
            }
            INFOS[k] = List::create(_["tag"] = IFD_tag,
                                    _["typ"] = IFD_type, 
                                    _["siz"] = IFD_count, 
                                    _["val"] = IFD_value, 
                                    _["byt"] = IFD_bytes, 
                                    _["len"] = IFD_count * multi[IFD_type], 
                                    _["off"] = IFD_off,
                                    _["map"] = IFD_map);
          }
        }
      }
      
      char buf_next [4];
      unsigned int next;
      fi.seekg(pos, std::ios::beg);
      fi.read((char*)&buf_next, sizeof(buf_next));
      memcpy(&next, buf_next, sizeof(next));
      if(swap) next = bytes_swap(next); 
      INFOS.names() = NAMES;
      List out = List::create(_["tags"] = INFOS,
                              _["infos"] = List::create(
                                _["IMAGE_LENGTH"] = IMAGE_LENGTH,
                                _["IMAGE_WIDTH"] = IMAGE_WIDTH,
                                // _["TILE_LENGTH"] = TILE_LENGTH,
                                // _["TILE_WIDTH"] = TILE_WIDTH,
                                _["OBJECT_ID"] = OBJECT_ID,
                                // _["OBJECT_ID_2"] = OBJECT_ID_2,
                                _["COMPRESSION"] = COMPRESSION,
                                _["TYPE"] = TYPE,
                                _["STRIP_OFFSETS"] = STRIP_OFFSETS,
                                _["STRIP_BYTE_COUNTS"] = STRIP_BYTE_COUNTS,
                                // _["TILE_OFFSETS"] = TILE_OFFSETS,
                                // _["TILE_BYTE_COUNTS"] = TILE_BYTE_COUNTS,
                                _["BG_MEAN"] = BG_MEAN,
                                _["BG_STD"] = BG_STD),
                                _["curr_IFD_offset"] = offset,
                                _["next_IFD_offset"] = next);
      out.attr("class") = "IFC_ifd";
      fi.close();
      return out;
    }
    catch(std::exception &ex) {
      fi.close();
      forward_exception_to_r(ex);
    }
    catch(...) { 
      Rcpp::stop("cpp_getTAGS: c++ exception (unknown reason)"); 
    }
  }
  else {
    Rcerr << "cpp_getTAGS: Unable to open " << fname << std::endl;
    Rcpp::stop("cpp_getTAGS: Unable to open file");
  }
  return List::create(_["tags"] = NA_REAL,
                      _["infos"] = List::create(
                        _["IMAGE_LENGTH"] = NA_REAL,
                        _["IMAGE_WIDTH"] = NA_REAL,
                        // _["TILE_LENGTH"] = NA_REAL,
                        // _["TILE_WIDTH"] = NA_REAL,
                        _["OBJECT_ID"] = NA_REAL,
                        // _["OBJECT_ID_2"] = NA_REAL,
                        _["COMPRESSION"] = NA_REAL,
                        _["TYPE"] = NA_REAL,
                        _["STRIP_OFFSETS"] = NA_REAL,
                        _["STRIP_BYTE_COUNTS"] = NA_REAL,
                        // _["TILE_OFFSETS"] = NA_REAL,
                        // _["TILE_BYTE_COUNTS"] = NA_REAL,
                        _["BG_MEAN"] = NA_REAL,
                        _["BG_STD"] = NA_REAL),
                        _["curr_IFD_offset"] = NA_REAL,
                        _["next_IFD_offset"] = NA_REAL);
}

//' @title IFC_offsets Computation with Object Identification
//' @name cpp_getoffsets_wid
//' @description
//' Returns offsets of the IFD (Image Field Directory) within a TIFF file.
//' @param fname string, path to file.
//' @param obj_count R_len_t, numbers of objects present in the file. Default is 0.
//' If obj_count <= 0 then progress_bar is forced to false.
//' @param display_progress bool, whether to display a progress bar. Default is false.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @source TIFF 6.0 specifications available at \url{https://www.adobe.io/open/standards/TIFF.html}
//' @return a list of integer vectors with OBJECT_ID, TYPE and OFFSET of IFDs found.
//' @keywords internal
////' @export
// [[Rcpp::export]]
List cpp_getoffsets_wid(std::string fname, R_len_t obj_count = 0, bool display_progress = false, bool verbose = false) {
  bool swap = false;
  std::string endianness = cpp_checkTIFF(fname);
  if(endianness == "big") swap = true;
  
  std::ifstream fi(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
  if (fi.is_open()) {
    try{
      fi.seekg(0, std::ios::end);
      if(obj_count == 0) display_progress = false;
      char buf_offset [4];
      uint32_t offset = 4;
      
      if(verbose) Rcout << "Extracting offsets from " << fname << std::endl;
      Progress p(obj_count * 2 + 1, display_progress);
      fi.seekg(offset, std::ios::beg);
      fi.read((char*)&buf_offset, sizeof(buf_offset));
      memcpy(&offset, buf_offset, sizeof(offset));
      if(swap) offset = bytes_swap(offset);
      if(!offset) {
        Rcerr << "cpp_getoffsets_wid: No IFD offsets found in\n" << fname << std::endl;
        Rcpp::stop("cpp_getoffsets_wid: No IFD offsets found");
      }
      
      IntegerVector out_obj, out_typ, out_off;
      IntegerVector obj = IntegerVector::create(NA_INTEGER);
      IntegerVector typ = IntegerVector::create(NA_INTEGER);

      uint32_t trunc_bytes = 8;
      bool force_trunc = true;
      while(offset){
        p.increment();
        List IFD = cpp_getTAGS(fname, offset, verbose, trunc_bytes, force_trunc);
        offset = as<uint32_t>(IFD["next_IFD_offset"]);
        List infos = IFD["infos"];
        
        if(iNotisNULL(infos["OBJECT_ID"])) {
          obj[0] = as<int32_t>(infos["OBJECT_ID"]); // min 0?, max ?
        }
        if(iNotisNULL(infos["TYPE"])) {
          typ[0] = as<int32_t>(infos["TYPE"]); // [0,3]
        }
        out_obj.push_back(obj[0]);
        out_typ.push_back(typ[0]);
        out_off.push_back(as<uint32_t>(IFD["curr_IFD_offset"]));
        
        if(Progress::check_abort()) {
          Rcerr << "cpp_getoffsets_wid: Interrupted by user" << std::endl;
          Rcpp::stop("cpp_getoffsets_wid: Interrupted by user");
        }
      }
      List out = List::create(_["OBJECT_ID"] = out_obj,
                              _["TYPE"] = out_typ,
                              _["OFFSET"] = out_off);
      fi.close();
      return out;
    }
    catch(std::exception &ex) {	
      fi.close();
      forward_exception_to_r(ex);
    }
    catch(...) { 
      Rcpp::stop("cpp_getoffsets_wid: c++ exception (unknown reason)"); 
    }
  }
  else {
    Rcerr << "cpp_getoffsets_wid: Unable to open " << fname << std::endl;
    Rcpp::stop("cpp_getoffsets_wid: Unable to open file");
  }
  return List::create(List::create(_["OBJECT_ID"] = NA_INTEGER,
                                   _["TYPE"] = NA_INTEGER,
                                   _["OFFSET"] = NA_INTEGER));
}

//' @title Checksum for RIF/CIF
//' @name cpp_checksum
//' @description
//' Computes sum of img IFDs (Image Field Directory) offsets of objects 0, 1, 2, 3 and 4.
//' @param fname string, path to file.
//' @source TIFF 6.0 specifications available at \url{https://www.adobe.io/open/standards/TIFF.html}
//' @return an integer vector with offsets of IFDs found.
//' @keywords internal
////' @export
// [[Rcpp::export]]
size_t cpp_checksum(std::string fname) {
  bool swap = false;
  std::string endianness = cpp_checkTIFF(fname);
  if(endianness == "big") swap = true;
  
  std::ifstream fi(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
  if (fi.is_open()) {
    try {
      fi.seekg(0, std::ios::end);
      char buf_offset [4];
      uint32_t offset = 4;
      IntegerVector obj = IntegerVector::create(0,1,2,3,4);
      IntegerVector found;
      uint8_t count = 0;
      size_t out = 0;
      
      fi.seekg(offset, std::ios::beg);
      fi.read((char*)&buf_offset, sizeof(buf_offset));
      memcpy(&offset, buf_offset, sizeof(offset));
      if(swap) offset = bytes_swap(offset);
      if(!offset) {
        Rcerr << "cpp_checksum: No IFD offsets found in\n" << fname << std::endl;
        Rcpp::stop("cpp_checksum: No IFD offsets found");
      }
      while(offset && (count < 5)){
        List IFD = cpp_getTAGS(fname, offset, false, 8, true);
        List infos = IFD["infos"];
        offset = as<uint32_t>(IFD["next_IFD_offset"]);
        if(iNotisNULL(infos["OBJECT_ID"])) {
          int32_t id = as<int32_t>(infos["OBJECT_ID"]);
          if(is_true(any(obj == id)) && !is_true(any(found == id))) {
            found.push_back(id);
            count += 1; 
            out += as<uint32_t>(IFD["curr_IFD_offset"]);
          } else {
            Rcpp::warning("cpp_checksum: raw object are not stored in expected order");
          }
        }
      }
      fi.close();
      return out;
    }
    catch(std::exception &ex) {	
      fi.close();
      forward_exception_to_r(ex);
    }
    catch(...) { 
      Rcpp::stop("cpp_checksum: c++ exception (unknown reason)"); 
    }
  }
  else {
    Rcerr << "cpp_checksum: Unable to open " << fname << std::endl;
    Rcpp::stop("cpp_checksum: Unable to open file");
  }
  return 0;
}

//' @title RLE Decompression
//' @name cpp_rle_Decomp
//' @description
//' Operates RLE decompression of compressed image stored in TIFF file.
//' @param fname string, path to file.
//' @param offset uint32_t, position of the beginning of compressed image.
//' @param nbytes uint32_t, number of bytes of compressed image.
//' @param imgWidth R_len_t, Width of the decompressed image. Default is 1.
//' @param imgHeight R_len_t, Height of the decompressed image. Default is 1.
//' @param nb_channels R_len_t, number of channels of the decompressed image. Default is 1.
//' @param removal uint8_t, object removal method. Default is no removal. Otherwise, if\cr
//' -1, for clipped removal: height OR width clipped pixels will be set to -1.\cr
//' -2, height clipped removal: height clipped pixels will be set to -1.\cr
//' -3, width clipped removal: width clipped pixels will be set to -1.\cr
//' -4, only keep background: background pixels will be set to 1 and all others to -1.\cr
//' -5, only keep foreground: foreground pixels will be set to 1 and all others to -1
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @details
//' BSD implementations of Bio-Formats readers and writers
//' %%
//' Copyright (C) 2005 - 2017 Open Microscopy Environment:
//'   - Board of Regents of the University of Wisconsin-Madison
//'   - Glencoe Software, Inc.
//'   - University of Dundee
//' %%
//' Redistribution and use in source and binary forms, with or without
//' modification, are permitted provided that the following conditions are met:
//' 
//' 1. Redistributions of source code must retain the above copyright notice,
//'    this list of conditions and the following disclaimer.
//' 2. Redistributions in binary form must reproduce the above copyright notice,
//'     this list of conditions and the following disclaimer in the documentation
//'     and/or other materials provided with the distribution.
//'  
//' THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//' IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//' ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
//' LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//' CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//' SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//' INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//' CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//' ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//' POSSIBILITY OF SUCH DAMAGE.
//' @source For image decompression, Lee Kamentsky's code porting from \url{https://github.com/openmicroscopy/bioformats/blob/4146b9a1797501f0fec7d6cfe69124959bff96ee/components/formats-bsd/src/loci/formats/in/FlowSightReader.java}
//' cited in \url{http://linkinghub.elsevier.com/retrieve/pii/S1046-2023(16)30291-2}
//' @keywords internal
////' @export
// [[Rcpp::export]]
List cpp_rle_Decomp (std::string fname, 
                     const uint32_t offset,
                     const uint32_t nbytes,
                     const R_len_t imgWidth = 1,
                     const R_len_t imgHeight = 1,
                     const R_len_t nb_channels = 1,
                     const uint8_t removal = 0,
                     const bool verbose = false) {
  if(nb_channels * imgWidth * imgHeight) {
    List out(nb_channels);
    R_len_t tile_width = imgWidth / nb_channels;
    std::ifstream fi(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
    if (fi.is_open()) {
      try{
        fi.seekg(0, std::ios::end);
        unsigned int filesize = fi.tellg();
        if(offset > (filesize - nbytes)) {
          Rcerr << "cpp_rle_Decomp: @offset:" << offset << " points to outside of\n" << fname << std::endl;
          Rcpp::stop("cpp_rle_Decomp: RLE image offset is higher than file size");
        }
        if(verbose) {
          Rcout << fname << std::endl;
          Rcout << "Extracting " << nbytes << " Bytes BitMask image @offset:" << offset << std::endl;
        }
        fi.seekg(offset, std::ios::beg);
        std::vector<char> buf_image(nbytes);
        fi.read(buf_image.data(), nbytes);
        
        IntegerMatrix img(imgWidth,imgHeight);
        uint32_t off, j, k, L = imgWidth * imgHeight, runLength = 0;
        
        int value;
        switch(removal) {
        case 1: { // clipped removal
          for(k = 0; k < nbytes; k++) {
          value = buf_image[k++];
          if(value > 1) value = -1;
          off = runLength;
          runLength = off + (buf_image[k] & 0xff) + 1;
          if (runLength > L) {
            Rcerr << "cpp_rle_Decomp: Buffer overrun in\n" << fname << std::endl;
            Rcpp::stop("cpp_rle_Decomp: Buffer overrun");
          }
          for(j = off; j < runLength; j++) img[j] = value;
        }
          break;
        }
        case 2: { // height clipped removal
          for(k = 0; k < nbytes; k++) {
          value = buf_image[k++];
          if(value == 2) value = -1;
          off = runLength;
          runLength = off + (buf_image[k] & 0xff) + 1;
          if (runLength > L) {
            Rcerr << "cpp_rle_Decomp: Buffer overrun in\n" << fname << std::endl;
            Rcpp::stop("cpp_rle_Decomp: Buffer overrun");
          }
          for(j = off; j < runLength; j++) img[j] = value;
        }
          break;
        }
        case 3: { // width clipped removal
          for(k = 0; k < nbytes; k++) {
          value = buf_image[k++];
          if(value == 3) value = -1;
          off = runLength;
          runLength = off + (buf_image[k] & 0xff) + 1;
          if (runLength > L) {
            Rcerr << "cpp_rle_Decomp: Buffer overrun in\n" << fname << std::endl;
            Rcpp::stop("cpp_rle_Decomp: Buffer overrun");
          }
          for(j = off; j < runLength; j++) img[j] = value;
        }
          break;
        }
        case 4: { // only keep background
          for(k = 0; k < nbytes; k++) {
          value = buf_image[k++];
          value = (value == 0) ? 1:-1;
          off = runLength;
          runLength = off + (buf_image[k] & 0xff) + 1;
          if (runLength > L) {
            Rcerr << "cpp_rle_Decomp: Buffer overrun in\n" << fname << std::endl;
            Rcpp::stop("cpp_rle_Decomp: Buffer overrun");
          }
          for(j = off; j < runLength; j++) img[j] = value;
        }
          break;
        }
        case 5: { // only keep non clipped foreground
          for(k = 0; k < nbytes; k++) {
          value = buf_image[k++];
          if(value != 1) value = -1;
          off = runLength;
          runLength = off + (buf_image[k] & 0xff) + 1;
          if (runLength > L) {
            Rcerr << "cpp_rle_Decomp: Buffer overrun in\n" << fname << std::endl;
            Rcpp::stop("cpp_rle_Decomp: Buffer overrun");
          }
          for(j = off; j < runLength; j++) img[j] = value;
        }
          break;
        }
        default: { // no removal 
          for(k = 0; k < nbytes; k++) {
          value = buf_image[k++];
          off = runLength;
          runLength = off + (buf_image[k] & 0xff) + 1;
          if (runLength > L) {
            Rcerr << "cpp_rle_Decomp: Buffer overrun in\n" << fname << std::endl;
            Rcpp::stop("cpp_rle_Decomp: Buffer overrun");
          }
          for(j = off; j < runLength; j++) img[j] = value;
        }
          break;
        }
        }
        fi.close();
        IntegerMatrix timg = transpose(img);
        for(R_len_t i = 0; i < nb_channels; i++) {
          out[i] = timg(_, Rcpp::Range(tile_width * i, tile_width * (i+1) - 1));
        }
        return out;
      }
      catch(std::exception &ex) {	
        fi.close();
        forward_exception_to_r(ex);
      }
      catch(...) { 
        Rcpp::stop("cpp_rle_Decomp: c++ exception (unknown reason)"); 
      }
    }
    else {
      Rcerr << "cpp_rle_Decomp: Unable to open " << fname << std::endl;
      Rcpp::stop("cpp_rle_Decomp: Unable to open file");
    }
  } else {
    Rcerr << "cpp_rle_Decomp: imgWidth, imgHeight and nb_channels should be >0" << std::endl;
    Rcpp::stop("cpp_rle_Decomp: imgWidth, imgHeight and nb_channels should be >0");    
  }
  return R_NilValue;
}

//' @title GRAY Decompression
//' @name cpp_gray_Decomp
//' @description
//' Operates GrayScale decompression of compressed image stored in TIFF file.
//' @param fname string, path to file.
//' @param offset uint32_t, position of the beginning of compressed image.
//' @param nbytes uint32_t, number of bytes of compressed image.
//' @param imgWidth R_len_t, Width of the decompressed image. Default is 1.
//' @param imgHeight R_len_t, Height of the decompressed image. Default is 1.
//' @param nb_channels R_len_t, number of channels of the decompressed image. Default is 1.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @details
//' BSD implementations of Bio-Formats readers and writers
//' %%
//' Copyright (C) 2005 - 2017 Open Microscopy Environment:
//'   - Board of Regents of the University of Wisconsin-Madison
//'   - Glencoe Software, Inc.
//'   - University of Dundee
//' %%
//' Redistribution and use in source and binary forms, with or without
//' modification, are permitted provided that the following conditions are met:
//' 
//' 1. Redistributions of source code must retain the above copyright notice,
//'    this list of conditions and the following disclaimer.
//' 2. Redistributions in binary form must reproduce the above copyright notice,
//'    this list of conditions and the following disclaimer in the documentation
//'    and/or other materials provided with the distribution.
//' 
//' THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//' IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//' ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
//' LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//' CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//' SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//' INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//' CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//' ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//' POSSIBILITY OF SUCH DAMAGE.
//' @source For image decompression, Lee Kamentsky's code porting from \url{http://github.com/openmicroscopy/bioformats/blob/4146b9a1797501f0fec7d6cfe69124959bff96ee/components/formats-bsd/src/loci/formats/in/FlowSightReader.java}\cr
//' cited in \url{http://linkinghub.elsevier.com/retrieve/pii/S1046-2023(16)30291-2}
//' @keywords internal
////' @export
// [[Rcpp::export]]
List cpp_gray_Decomp (std::string fname, 
                      const uint32_t offset, 
                      const uint32_t nbytes, // force number of bytes to be less than 65,535
                      const R_len_t imgWidth = 1, 
                      const R_len_t imgHeight = 1, 
                      const R_len_t nb_channels = 1,
                      const bool verbose = false) {
  if(nb_channels * imgWidth * imgHeight) {
    List out(nb_channels);
    R_len_t tile_width = imgWidth / nb_channels;
    std::ifstream fi(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
    if(fi.is_open()) {
      try {
        fi.seekg(0, std::ios::end);
        size_t filesize = fi.tellg();
        if(offset > (filesize - nbytes)) {
          Rcerr <<  "cpp_gray_Decomp: @offset:" << offset << " points to outside of\n" << fname  << std::endl;
          Rcpp::stop("cpp_gray_Decomp: GrayScale image offset is higher than file size");
        }
        if(verbose) {
          Rcout << fname << std::endl;
          Rcout << "Extracting " << nbytes << " Bytes GreyScale image @offset:" << offset << std::endl;
        }
        fi.seekg(offset, std::ios::beg);
        std::vector<char> buf_image(nbytes);
        fi.read(buf_image.data(), nbytes);
        
        IntegerVector lastRow(imgWidth + 1);
        IntegerMatrix img(imgHeight, imgWidth + 1);
        bool odd = false;
        uint32_t k = 0;
        int shift, value, nibble;
        
        for(R_len_t y = 0 ; y < imgHeight ; y++) {
          for(R_len_t x = 1 ; x <= imgWidth ; x++) {
            shift = 0;
            value = 0;
            nibble = -1;
            while((nibble & 0x8)) {
              nibble = odd ? buf_image[k++] >> 4 : buf_image[k] & 0xf;
              odd = !odd;
              value += (nibble & 0x7) << shift;
              shift += 3;
            }
            if(nibble & 0x4) value |= - (1 << shift);
            lastRow[x] += value;
            img(y,x) = img(y,x - 1) + lastRow[x];
          }
        }
        fi.close();
        for(R_len_t i = 0; i < nb_channels; i++) {
          out[i] = img(_, Rcpp::Range(1 + tile_width * i, tile_width * (i+1)));
        }
        return out;
      }
      catch(std::exception &ex) {	
        fi.close();
        forward_exception_to_r(ex);
      }
      catch(...) { 
        Rcpp::stop("cpp_gray_Decomp: c++ exception (unknown reason)"); 
      }
    }
    else {
      Rcerr << "cpp_gray_Decomp: Unable to open " << fname << std::endl;
      Rcpp::stop("cpp_gray_Decomp: Unable to open file");
    }
  } else {
    Rcerr << "cpp_gray_Decomp: imgWidth, imgHeight and nb_channels should be >0" << std::endl;
    Rcpp::stop("cpp_gray_Decomp: imgWidth, imgHeight and nb_channels should be >0");    
  }
  return R_NilValue;
}


//' @title IFC_object Decompression
//' @name cpp_decomp
//' @description
//' Operates decompression of compressed image stored in TIFF file.
//' @param fname string, path to file.
//' @param offset uint32_t, position of the beginning of compressed image.
//' @param nbytes uint32_t, number of bytes of compressed image.
//' @param imgWidth R_len_t, Width of the decompressed image. Default is 1.
//' @param imgHeight R_len_t, Height of the decompressed image. Default is 1.
//' @param nb_channels R_len_t, number of channels of the decompressed image. Default is 1.
//' @param removal uint8_t, object removal method. Only apply for 30818 compression. Default is 0.\cr
//' -1, for clipped removal: height OR width clipped pixels will be set to -1.\cr
//' -2, height clipped removal: height clipped pixels will be set to -1.\cr
//' -3, width clipped removal: width clipped pixels will be set to -1.\cr
//' -4, only keep background: background pixels will be set to 1 and all others to 0.\cr
//' -5, only keep foreground: foreground pixels will be set to 1 and all others to 0.
//' @param compression uint32_t, compression algorithm used. Default is 30818.
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @details
//' BSD implementations of Bio-Formats readers and writers
//' %%
//' Copyright (C) 2005 - 2017 Open Microscopy Environment:
//'   - Board of Regents of the University of Wisconsin-Madison
//'   - Glencoe Software, Inc.
//'   - University of Dundee
//' %%
//' Redistribution and use in source and binary forms, with or without
//' modification, are permitted provided that the following conditions are met:
//' 
//' 1. Redistributions of source code must retain the above copyright notice,
//'    this list of conditions and the following disclaimer.
//' 2. Redistributions in binary form must reproduce the above copyright notice,
//'    this list of conditions and the following disclaimer in the documentation
//'    and/or other materials provided with the distribution.
//' 
//' THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//' IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//' ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
//' LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//' CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//' SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//' INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//' CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//' ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//' POSSIBILITY OF SUCH DAMAGE.
//' @source For image decompression, Lee Kamentsky's code porting from \url{http://github.com/openmicroscopy/bioformats/blob/4146b9a1797501f0fec7d6cfe69124959bff96ee/components/formats-bsd/src/loci/formats/in/FlowSightReader.java}\cr
//' cited in \url{http://linkinghub.elsevier.com/retrieve/pii/S1046-2023(16)30291-2}
//' @keywords internal
////' @export
// [[Rcpp::export]]
List cpp_decomp (std::string fname, 
                 const uint32_t offset, 
                 const uint32_t nbytes, 
                 const R_len_t imgWidth = 1, 
                 const R_len_t imgHeight = 1, 
                 const R_len_t nb_channels = 1,
                 const uint8_t removal = 0,
                 const uint32_t compression = 1,
                 const bool verbose = false) {
  switch(compression) {
  case 30817: return cpp_gray_Decomp(fname, offset, nbytes, imgWidth, imgHeight, nb_channels, verbose);
  case 30818: return cpp_rle_Decomp(fname, offset, nbytes, imgWidth, imgHeight, nb_channels, removal, verbose);
  }
  Rcerr << "cpp_decomp: can't deal with compression format:" << compression << std::endl;
  Rcpp::stop("cpp_decomp: can't deal with compression format");   
  return R_NilValue;
}

//' @title Matrix Normalization
//' @name cpp_normalize
//' @description
//' Normalizes a finite matrix to [0,1]
//' @param mat a finite NumericMatrix.
//' @param input_range a finite NumericVector, sets the range of the input intensity values. Default is c(0,4095).
//' values outside this range are clipped.
//' @param full_range if full_range is TRUE, then input_range will be set to 'c(0,4095)' and gamma forced to 1. Default is false.
//' @param force_range if force_range is TRUE, then input_range will be adjusted to mat range and gamma forced to 1. Default is false.\cr
//' Note that this parameter takes the precedence over 'input_range' and 'full_range'.
//' @param gamma correction. Default is 1, for no correction.
//' @keywords internal
////' @export
// [[Rcpp::export]]
NumericMatrix cpp_normalize (const NumericMatrix mat, 
                             const NumericVector input_range = NumericVector::create(0.0,4095.0),
                             const bool full_range = false,
                             const bool force_range = false, 
                             const double gamma = 1.0) {
  NumericMatrix xx = no_init_matrix(mat.nrow(), mat.ncol());
  NumericVector ran(2);
  double gam = gamma;
  if(force_range) {
    gam = 1.0;
    ran = cpp_check_range(cpp_check_range(mat)); // twice to ensure that mat is at least of length 2
  } else {
    if(full_range) {
      gam = 1.0;
      ran[0] = 0.0;
      ran[1] = 4095.0;
    } else {
      ran = cpp_check_range(cpp_check_range(input_range)); // twice to ensure that input_range is at least of length 2
    }
  }

  double diff = ran[1] - ran[0];
  if(gamma == 1) {
    for(R_len_t i = 0; i < mat.size(); i++) {
      if(mat[i] <= ran[0]) {
        xx[i] = 0.0;
        continue;
      }
      if(mat[i] >= ran[1]) {
        xx[i] = 1.0;
        continue;
      }
      xx[i] = (mat[i] - ran[0])/diff;
    }
  } else {
    for(R_len_t i = 0; i < mat.size(); i++) {
      if(mat[i] <= ran[0]) {
        xx[i] = 0.0;
        continue;
      }
      if(mat[i] >= ran[1]) {
        xx[i] = 1.0;
        continue;
      }
      xx[i] = pow((mat[i] - ran[0])/diff, gam);
    }
  }
  return xx;
}

//' @title Matrix Cleanser
//' @name cpp_cleanse
//' @description
//' Replaces values in matrix mat according to mask msk.
//' Depending of add_noise parameter values will be replaced with noise or not.
//' @param mat a NumericMatrix.
//' @param msk a LogicalMatrix.
//' @param add_noise logical, if true adds normal noise.
//' Rcpp::Rf_rnorm(bg, sd) function is used. Default is true.
//' @param bg double, mean value of the background added if add_noise is true. Default is 0.
//' @param sd double, standard deviation of the background added if add_noise is true. Default is 0.
//' @return a NumericMatrix with replaced according to msk.
//' @keywords internal
////' @export
// [[Rcpp::export]]
NumericMatrix cpp_cleanse (const NumericMatrix mat, 
                           const LogicalMatrix msk, 
                           const bool add_noise = true, 
                           const double bg = 0.0, const double sd = 0.0) {
  if(!(msk.ncol() == mat.ncol()) && (msk.nrow() == mat.nrow())) Rcpp::stop("cpp_cleanse: mat and msk should have same dimensions");
  NumericMatrix out = no_init_matrix(mat.nrow(), mat.ncol());
  R_len_t i = 0;
  if(add_noise) {
    for(; i < out.size(); i++) out[i] = msk[i] ? Rf_rnorm(bg, sd):mat[i];
  } else {
    for(; i < out.size(); i++) out[i] = msk[i] ? bg:mat[i];
  }
  return out;
}

//' @title Equal Sized Matrix to Matrix Writer According to Mask
//' @name cpp_mask
//' @description
//' Writes matrix B in matrix A according to mask.
//' If mask is not 0 B is written, A otherwise.
//' @param A a NumericMatrix.
//' @param B a NumericMatrix.
//' @param mask a NumericMatrix.
//' @keywords internal
////' @export
// [[Rcpp::export]]
NumericMatrix cpp_mask (const NumericMatrix A,
                        const NumericMatrix B,
                        const NumericMatrix mask) {
  R_len_t ar = A.nrow(), ac = A.ncol();
  if((B.ncol() != ac) || (mask.ncol() != ac) || (B.nrow() != ar) || (mask.nrow() != ar) ) Rcpp::stop("cpp_mask: A, B and mask should have same dimensions");
  Rcpp::NumericMatrix out(ar, ac);
  for(R_len_t i = 0; i < A.size(); i++) out[i] = mask[i] ? B[i] : A[i];
  return out;
}

//' @title Matrix to Matrix Writer According to Mask with Offsets
//' @name cpp_mark
//' @description
//' Writes matrix B in matrix A according to mask.
//' @param A a NumericMatrix.
//' @param B a NumericMatrix.
//' @param mask a NumericMatrix.
//' @param xoff x offset in A to start writing B.
//' @param yoff x offset in A to start writing B.
//' @param invert a logical. Default is false.
//' When false, the default, values of B are written into A when mask is not 0.
//' When true, values of 1-B are written into A when mask is not 0.
//' @keywords internal
////' @export
// [[Rcpp::export]]
NumericMatrix cpp_mark (const NumericMatrix A,
                        const NumericMatrix B,
                        const NumericMatrix mask,
                        const R_len_t xoff = 0,
                        const R_len_t yoff = 0,
                        const bool invert = false) {
  R_len_t bc = B.ncol();
  R_len_t br = B.nrow();
  if((A.ncol() < (bc + xoff)) || (A.nrow() < (br + yoff))) Rcpp::stop("cpp_mark: A should be at least of same dimensions as B + offsets");
  if((mask.ncol() < bc) || (mask.nrow() < br)) Rcpp::stop("cpp_mark: mask should be at least of same dimensions as B");
  Rcpp::NumericMatrix out = Rcpp::clone(A);
  if(invert) {
    for(R_len_t y = 0; y < br; y++) for(R_len_t x = 0; x < bc; x++) if(mask(y,x)) out(y+yoff,x+xoff) = std::fabs(1-B(y,x));
  } else {
    for(R_len_t y = 0; y < br; y++) for(R_len_t x = 0; x < bc; x++) if(mask(y,x)) out(y+yoff,x+xoff) = B(y,x);
  }
  return out;
}

//' @title Matrix Transformation
//' @name cpp_transform
//' @description
//' Function to normalize, colorize and add background to images.
//' @param mat NumericMatrix.
//' @param color NumericVector, whose members are h,s,v color.
//' This vector has to be named with 1st name being the name of this color.
//' @param msk LogicalMatrix.
//' @param size a length 2 IntegerVector, of final dimensions (height,width) of the image. Default is 0,0 for no change.
//' @param mode string, color mode export. Either "rgb", "gray" or "raw". Default is "raw".
//' @param type uint16_t image object type.
//' Rcpp::Rf_rnorm(bg, sd) function is used. Default is true.
//' @param cleanse bool, wheteher to cleanse image or not. Cleanse will only apply when set to true and 'msk' is not R_NilValue
//' @param input_range a finite NumericVector, only apply when mode is not "raw", sets the range of the input intensity values. Default is c(0,4095).
//' values exceeding this range are clipped.
//' @param add_noise bool, if true adds normal noise.
//' @param bg double, mean value of the background added if add_noise is true. Default is 0.
//' @param sd double, standard deviation of the background added if add_noise is true. Default is 0.
//' @param full_range bool, only apply when mode is not "raw", if full_range is TRUE, then input_range will be set to 'c(0,4095)' and gamma forced to 1. Default is false.
//' @param force_range bool, only apply when mode is not "raw", if force_range is TRUE, then input_range will be adjusted to mat range and gamma forced to 1. Default is false.\cr
//' Note that this parameter takes the precedence over 'input_range' and 'full_range'.
//' @param gamma correction. Default is 1, for no correction.
//' @keywords internal
////' @export
// [[Rcpp::export]]
NumericVector cpp_transform(const NumericMatrix mat,
                            const NumericVector color,
                            const LogicalMatrix msk,
                            const IntegerVector size = IntegerVector::create(0,0),
                            const std::string mode = "raw",
                            const uint16_t type = 2,
                            const bool cleanse = false,
                            const NumericVector input_range = NumericVector::create(0.0,4095.0),
                            const bool add_noise = true,
                            const double bg = 0.0,
                            const double sd = 0.0,
                            const bool full_range = false,
                            const bool force_range = false,
                            const double gamma = 1.0) {
  NumericMatrix foo;
  if(cleanse) {
    foo = cpp_cleanse(mat, msk, add_noise, bg, sd);
  } else {
    foo = Rcpp::clone(mat);
  }
  foo = cpp_resize(foo, size[0], size[1], add_noise, bg, sd);
  Rcpp::CharacterVector col_name = wrap(color.attr("names"));
  if(mode != "raw") {
    foo = cpp_normalize(foo, input_range, full_range, force_range, gamma);
  }
  if(mode == "rgb") {
    NumericVector bar = cpp_M_HSV2RGB(foo, color[0], color[1]);
    bar.attr("dim") = Rcpp::Dimension(foo.nrow(), foo.ncol(), 3);
    bar.attr("input_range") = input_range;
    bar.attr("full_range") = full_range;
    bar.attr("force_range") = force_range;
    bar.attr("gamma") = gamma;
    bar.attr("color") = col_name[0];
    bar.attr("mode") = mode;
    bar.attr("RAW") = mat;
    bar.attr("BG_MEAN") = (bg >= 0.0) ? bg:0.0; // negative values (i.e. -1) are used for removal of non masked objects
    bar.attr("BG_STD") = sd;
    if(type == 2) {
      bar.attr("class") = "IFC_img";
    } else {
      bar.attr("class") = "IFC_msk";
    }
    return bar;
  }
  foo.attr("input_range") = input_range;
  foo.attr("full_range") = full_range;
  foo.attr("force_range") = force_range;
  foo.attr("gamma") = gamma;
  foo.attr("color") = "Gray";
  foo.attr("mode") = mode;
  foo.attr("RAW") = mat;
  foo.attr("BG_MEAN") = (bg >= 0.0) ? bg:0.0; // negative values (i.e. -1) are used for removal of non masked objects
  foo.attr("BG_STD") = sd;
  if(type == 2) {
    foo.attr("class") = "IFC_img";
  } else {
    foo.attr("class") = "IFC_msk";
  }
  return foo;
}

//' @title IFC_object Extraction
//' @name cpp_extract
//' @description
//' Extracts object from ifd
//' @param fname string, path to file
//' @param ifd List, ifd information of class IFC_ifd
//' @param colors List of colors to use.
//' @param channels DataFrame, channels information
//' @param chan_to_extract IntegerVector, channels to extract
//' @param extract_msk uint8_t, type of masked to extract.\cr
//' - 0: no mask\cr
//' - 1: at least one clipped\cr
//' - 2: at least one masked\cr
//' - 3: at least one MC
//' @param mode string, color mode export. Either "rgb", "gray" or "raw". Default is "raw".
//' @param size a length 2 IntegerVector of final dimensions (height,width) of the image. Default is 0,0 for no change.\cr
//' @param verbose bool, whether to display information (use for debugging purpose). Default is false.
//' @keywords internal
////' @export
// [[Rcpp::export]]
List cpp_extract (std::string fname,
                  List ifd,
                  const List colors,
                  const DataFrame channels,
                  const IntegerVector chan_to_extract,
                  const uint8_t extract_msk = 0,
                  const std::string mode = "raw",
                  const IntegerVector size = IntegerVector::create(0,0),
                  bool verbose = false) {
  
  R_len_t nb_channels = channels.nrows();
  List infos = ifd["infos"];
  R_len_t iml = infos["IMAGE_LENGTH"];
  R_len_t imw = infos["IMAGE_WIDTH"];
  uint16_t typ = infos["TYPE"];
  uint32_t off = infos["STRIP_OFFSETS"];
  uint32_t byt = infos["STRIP_BYTE_COUNTS"];
  uint32_t com = infos["COMPRESSION"];
  
  CharacterVector _physicalChannel = channels["physicalChannel"];
  NumericVector _xmin = channels["xmin"];
  NumericVector _xmax = channels["xmax"];
  IntegerVector _removal = channels["removal"];
  LogicalVector _add_noise = channels["add_noise"];
  NumericVector _bg;
  NumericVector _sd;
  LogicalVector _full_range = channels["full_range"];
  LogicalVector _force_range = channels["force_range"];
  NumericVector _gamma = channels["gamma"];
  
  List out(chan_to_extract.size());
  switch(typ) {
  case 1: { // 1st offset is detected, nothing to extract
    for(R_len_t i = 0; i < chan_to_extract.size(); i++) {
    out[i] = NumericMatrix(iml, imw);
  }
    return out;
    break;
  }
  case 2: { // an image is detected use background mean and sd
    _bg = infos["BG_MEAN"];
    _sd = infos["BG_STD"];
    break;
  }
  case 3: { // a mask is detected some parameters are forced
    _xmin = rep(0, nb_channels);
    for(R_len_t i = 0; i < nb_channels; i++) _xmax[i] = _removal[i] ? 1:3;
    _bg = rep(0, nb_channels);
    _sd = rep(0, nb_channels);
    _gamma = rep(1, nb_channels);
    _add_noise = rep(false, nb_channels);
    _force_range = rep(false, nb_channels);
    break;
  }
  default: { // not allowed type
    Rcerr <<  "cpp_extract: trying to extract a unknow object";
    Rcpp::stop("cpp_extract: trying to extract a unknow object");
  }
  }

  // extract image
  List img = cpp_decomp(fname, off, byt, 
                        imw, iml, nb_channels,
                        0, com, verbose);
  LogicalMatrix msk_rm(iml, imw / nb_channels);
  if(extract_msk > 0) {
    List masks;
    LogicalMatrix MC(iml, imw / nb_channels);
    MC.fill(true);
    if(typ == 2) {
      List msk_ifd = cpp_getTAGS(fname, ifd["next_IFD_offset"], verbose, 8, true)["infos"];
      masks = cpp_decomp(fname, msk_ifd["STRIP_OFFSETS"], msk_ifd["STRIP_BYTE_COUNTS"], 
                         imw, iml, nb_channels, 
                         0, msk_ifd["COMPRESSION"], verbose);
      
    } else {
      masks = clone(img);
    }
    if(extract_msk == 3) {
      for(R_len_t i = 0; i < masks.length(); i++) {
        NumericMatrix CUR_M = clone(Rcpp::as<Rcpp::NumericMatrix>(masks[i]));
        for(R_len_t i_row = 0; i_row < iml; i_row ++) {
          MC(i_row, _) = (CUR_M(i_row, _) != 1) & MC(i_row, _);
        }
      }
    }
    masks.attr("names") = _physicalChannel;
    
    // transform extracted image according to user's settings with mask removal
    for(R_len_t i = 0; i < chan_to_extract.size(); i++) {
      uint8_t chan_idx = chan_to_extract[i];
      std::string cur_chan = as<std::string>(_physicalChannel[chan_idx]);
      double msk_bm = _bg[chan_idx];
      switch(_removal[chan_idx]) {
      case 1: {
        NumericMatrix CUR_M = clone(Rcpp::as<Rcpp::NumericMatrix>(masks[cur_chan]));
        for(R_len_t i_row = 0; i_row < iml; i_row ++) {
          msk_rm(i_row, _) = CUR_M(i_row, _) > 1;
        }
        break;
      }
      case 2: {
        NumericMatrix CUR_M = clone(Rcpp::as<Rcpp::NumericMatrix>(masks[cur_chan]));
        for(R_len_t i_row = 0; i_row < iml; i_row ++) {
          msk_rm(i_row, _) = CUR_M(i_row, _) != 1;
        }
        if((typ == 2) && (!_add_noise[chan_idx])) msk_bm = -1;
        break;
      }
      case 3: {
        msk_rm = MC;
        if((typ == 2) && (!_add_noise[chan_idx])) msk_bm = -1;
        break;
      }
      }
      out[i] = cpp_transform(img[chan_idx],
                             colors[chan_idx],
                                   msk_rm,
                                   size,
                                   mode,
                                   typ,
                                   _removal[chan_idx],
                                           NumericVector::create(_xmin[chan_idx],_xmax[chan_idx]),
                                           _add_noise[chan_idx],
                                                     msk_bm,
                                                     _sd[chan_idx],
                                                        _full_range[chan_idx],
                                                                   _force_range[chan_idx],
                                                                               _gamma[chan_idx]);
    }
  } else {
    for(R_len_t i = 0; i < chan_to_extract.size(); i++) {  
      // transform extracted image according to user's settings without mask removal
      uint8_t chan_idx = chan_to_extract[i];
      out[i] = cpp_transform(img[chan_idx],
                             colors[chan_idx],
                                   msk_rm,
                                   size,
                                   mode,
                                   typ,
                                   false,
                                   NumericVector::create(_xmin[chan_idx],_xmax[chan_idx]),
                                   _add_noise[chan_idx],
                                             _bg[chan_idx],
                                                _sd[chan_idx],
                                                   _full_range[chan_idx],
                                                              _force_range[chan_idx],
                                                                          _gamma[chan_idx]);
    }
  }
  return out;
}
