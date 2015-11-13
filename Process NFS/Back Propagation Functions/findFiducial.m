function [ pos ] = findFiducial( search_region,sample_size)

%%% Left Bottom
[E_lb,yyy, zzz]=reSquare(search_region.lb.fields, search_region.lb.y, search_region.lb.z);
[I,J] = findsubmax(E_lb,20);
pos.left_b(1)=(J-1)*sample_size+min(min(yyy));
pos.left_b(2)=(I-1)*sample_size+min(min(zzz));

%%% Right Bottom
[E_rb,yyy, zzz]=reSquare(search_region.rb.fields, search_region.rb.y, search_region.rb.z);
[I,J] = findsubmax(E_rb,20);
pos.right_b(1)=(J-1)*sample_size+min(min(yyy)); %y 
pos.right_b(2)=(I-1)*sample_size+min(min(zzz)); %z

%%% Left Top
[E_lt,yyy, zzz]=reSquare(search_region.lt.fields, search_region.lt.y, search_region.lt.z);
[I,J] = findsubmax(E_lt,20);
pos.left_t(1)=(J-1)*sample_size+min(min(yyy));
pos.left_t(2)=(I-1)*sample_size+min(min(zzz));

%%% Right Top
[E_rt,yyy, zzz]=reSquare(search_region.rt.fields, search_region.rt.y, search_region.rt.z);
[I,J] = findsubmax(E_rt,20);
pos.right_t(1)=(J-1)*sample_size+min(min(yyy));
pos.right_t(2)=(I-1)*sample_size+min(min(zzz));
 

end
