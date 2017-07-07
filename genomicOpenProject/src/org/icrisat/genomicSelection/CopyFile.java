/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.icrisat.genomicSelection;

import java.io.*;

/**
 *
 * @author Mohan
 * Copying the file from source directory to result directory
 */
public class CopyFile{
        public static void copy(String source, String destination) throws IOException {
        File source_file = new File(source);
        File desti_file = new File(destination);
        FileInputStream fis = null;
        FileOutputStream fos = null;    
        byte[] buffer;
        int byte_read;
        try {
            /* first check that file exists or not. */
            if (!source_file.exists() || !source_file.isFile()) 
            {
                throw new ClassException("No source file found : "+source);
            }
	    /* check that the file is readable or not. */
            if (!source_file.canRead()) {
                throw new ClassException("Source file is unreadable: "+source);
            }
            /* If the destination exists, make sure it is a 
            writeable file and ask before overwriting it. 
            If the destination doesn't exist, make sure the 
            directory exists and is writeable.*/
            if (desti_file.exists()) {
                if (desti_file.isFile()) {
                    DataInputStream in = 
                    new DataInputStream(System.in);
                    if (!desti_file.canWrite()) {
                        throw new ClassException("Destination is unwriteable : "+destination);
                    }
                    
                    System.out.flush();
                    String response = "y";
                    if (!response.equals("Y") && 
                         !response.equals("y")) {
                        throw new ClassException("Wrong Input.");
                    }
                } 
		else {
                   throw new ClassException("Destination is not a normal file: " + destination);
                }
            } 
	    else {
                File parentdir = parent(desti_file);
                if (!parentdir.exists()) {
                    throw new ClassException("No Destination directory exist: " + destination);
                }
                if (!parentdir.canWrite()) {
                    throw new ClassException("Destination directory is unwriteable: " + destination);
                }
            }
            /* Now we have checked all the things so we can 
            copy the file now.*/
            fis = new FileInputStream(source_file);
            fos = new FileOutputStream(desti_file);
            buffer = new byte[1024];
            while (true) {
                byte_read = fis.read(buffer);
                if (byte_read == -1) {
                    break;
                }
                fos.write(buffer, 0, byte_read);
            }
        } 
        /* Finally close the stream. */ 
	finally {
              fis.close();
	      fos.close();              
        } }
    /* File.getParent() can return null when the file is 
    specified without a directory or
	is in the root directory. This method handles those cases.*/
    private static File parent(File f) {
        String dirname = f.getParent();
        if (dirname == null) {
            if (f.isAbsolute()) {
                return new File(File.separator);
            } 
	    else {
                return new File(System.getProperty("user.dir"));
            }
        }
        return new File(dirname);
    }
    
    public static void main(String[] args) {
        if (args.length != 2) {
            System.err.println("Wrong argument entered.  Try Again......");
        } 
	else {
            try {
                copy(args[0], args[1]);
            } 
	    catch (IOException ex) {
                System.err.println(ex.getMessage());
            }
        }
    }
}
class ClassException extends IOException {
    public ClassException(String msg) {
        super(msg);
    }
}
