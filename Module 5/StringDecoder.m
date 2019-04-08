in_message_string = 'potato';
in_message_binary = reshape(logical(dec2bin(in_message_string) - 48).',1,[]);
display(in_message_binary);


out_message_binary = in_message_binary;


out_message_string = char(bin2dec(char(reshape(out_message_binary,7,[]).' + 48))).';
display(out_message_string);
