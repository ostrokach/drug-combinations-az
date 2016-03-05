DELIMITER $$
 
CREATE FUNCTION gene_pair_essentiality(gene_essentiality_1 smallint, gene_essentiality_2 smallint) RETURNS smallint
  DETERMINISTIC
BEGIN
  DECLARE result smallint;
 
  IF gene_essentiality_1 = -1 and gene_essentiality_2 = -1 then
    set result = 1;
  ELSEIF (gene_essentiality_1 = -1 and gene_essentiality_2 = 0) or
       (gene_essentiality_1 = 0 and gene_essentiality_2 = -1) then
    set result = 2;
  elseif gene_essentiality_1 = 0 and gene_essentiality_2 = 0 then
    set result = 3;
  elseif (gene_essentiality_1 = 1 and gene_essentiality_2 = 0) or
       (gene_essentiality_1 = 0 and gene_essentiality_2 = 1) then
    set result = 4;
  elseif gene_essentiality_1 = 1 and gene_essentiality_2 = 1 then
    set result = 5;
  elseif (gene_essentiality_1 = -1 and gene_essentiality_2 = 1) or
       (gene_essentiality_1 = 1 and gene_essentiality_2 = -1) then
    set result = 6;
  END IF;
 
 RETURN (result);
END
